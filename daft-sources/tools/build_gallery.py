#!/usr/bin/env python3
"""Publish the book's numbered figures, previews, source bundle, and notebooks."""
import argparse
import hashlib
import html
import io
import json
from pathlib import Path
import re
import zipfile

import pypdfium2 as pdfium
from pypdf import PdfReader

PUBLIC = 'https://mathematical-tours.github.io'
GITHUB = 'https://github.com/mathematical-tours/mathematical-tours.github.io/blob/master/daft-sources/'
CHAPTERS = {
    1: 'Fourier analysis on finite groups',
    2: 'Applications of Fourier duality',
    3: 'The discrete Fourier transform',
    4: 'Applications of the discrete Fourier transform',
    5: 'Transforms related to the DFT',
    6: 'Fourier transforms over finite fields and rings',
    7: 'Linear representations of finite groups',
}
TITLES = {
    'cyclic-characters': 'Characters of a cyclic group',
    'fourier-decomposition': 'Fourier decomposition',
    'circular-convolution': 'Circular convolution',
    'random-walks': 'Random walks on a cyclic group',
    'boolean-cube': 'The Boolean cube',
    'walsh-matrix': 'The Walsh matrix',
    'geometric-reciprocity': 'Quadratic reciprocity on a lattice',
    'haar-projections': 'Successive Haar projections',
    'haar-cascade': 'The Haar cascade',
    'haar-examples': 'Discrete Haar transforms',
    'walsh-orderings': 'Natural and sequency order',
    'walsh-spectrum': 'The Walsh spectrum',
    'walsh-signal-approximation': 'Walsh approximation of a signal',
    'walsh-basis-two-dimensions': 'Two-dimensional Walsh functions',
    'walsh-image-approximation': 'Walsh approximation of an image',
    'intermediate-walsh-transform': 'Intermediate Walsh transforms',
    'elementary-butterfly': 'An FFT butterfly',
    'eight-point-fft': 'An eight-point FFT',
    'bit-reversal-matrices': 'Bit-reversal matrices',
    'bit-reversal-permutation': 'The bit-reversal permutation',
    'linear-and-circular-convolution': 'Linear and circular convolution',
    'two-dimensional-dft': 'The two-dimensional DFT',
    'two-dimensional-convolution': 'Two-dimensional convolution',
    'fourier-eigenspace-projections': 'Fourier eigenspace projections',
    'intermediate-fourier-transform': 'Intermediate Fourier transforms',
    'trigonometric-interpolation': 'Trigonometric interpolation',
    'chebyshev-polynomials': 'Chebyshev polynomials',
    'chebyshev-interpolation': 'Choosing interpolation nodes',
    'fractional-differentiation': 'Fractional differentiation',
    'fractional-fourier-matrices': 'Fractional Fourier matrices',
    'fractional-fourier-gaussian': 'Fractional transforms of a Gaussian',
    'fourier-eigenvectors': 'Structured Fourier eigenvectors',
    'fourier-eigenbasis-matrix': 'A Fourier eigenbasis',
    'orthogonalizing-translates': 'Orthogonalizing cyclic translates',
    'sampling-and-padding': 'Sampling and zero padding',
    'gaussian-signal-smoothing': 'Gaussian smoothing of a signal',
    'filter-responses': 'Impulse, frequency, and step responses',
    'gaussian-image-smoothing': 'Gaussian smoothing of an image',
    'image-denoising': 'Image denoising',
    'polygon-iteration-regimes': 'Four regimes of polygon filtering',
    'polygon-filtering': 'Filtering a polygon',
    'fourier-descriptors': 'Fourier shape descriptors',
    'heat-evolution': 'The spectral heat equation',
    'poisson-reference': 'Solving the Poisson equation',
    'poisson-modified-boundary': 'Changing the boundary values',
    'polynomial-convolution': 'Multiplying polynomials by convolution',
    'heat-time-stepping': 'Time stepping for the heat equation',
    'template-matching': 'Template matching by correlation',
    'fourier-image-rotation': 'Rotating an image with Fourier shears',
    'low-pass-filters': 'Designing low-pass filters',
    'centered-bsplines': 'Centered B-splines',
    'cardinal-spline-and-sinc': 'Cardinal splines and sinc',
    'cubic-spline-boundaries': 'Boundary conditions for cubic splines',
    'hartley-and-fourier': 'Hartley and Fourier spectra',
    'notch-poles-and-zeros': 'Poles and zeros of a notch filter',
    'notch-transfer-surface': 'The notch-filter transfer function',
    'notch-responses': 'Notch-filter responses',
    'rlc-circuit': 'A series RLC circuit',
    'z-transform-contours': 'Contours for the z-transform',
    'scaled-fourier-signal': 'Scaled-frequency Fourier transforms',
    'spectral-leakage': 'Spectral leakage',
    'spectral-zoom-and-adjustment': 'Zooming into a spectrum',
    'fractional-hartley': 'Fractional Hartley transforms',
    'generalized-hartley': 'Generalized Hartley transforms',
    'scaled-fourier-image': 'Scaled-frequency image transforms',
    'quadrature-filter-responses': 'Quadrature filter responses',
    'hamming-cube': 'Hamming distance on the cube',
    'finite-field-wavelets': 'Wavelets over a finite field',
}


def write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(content, bytes):
        path.write_bytes(content)
    else:
        path.write_text(content, encoding='utf-8')


def group(text, start):
    """Read a balanced TeX group, leaving escaped braces intact."""
    assert text[start] == '{'
    depth = 1
    i = start + 1
    while i < len(text):
        if text[i] == '\\':
            i += 2
            continue
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                return text[start + 1:i], i + 1
        i += 1
    raise ValueError('Unclosed TeX group')


def labels(root, path, seen=None):
    seen = set() if seen is None else seen
    if path in seen:
        return {}
    seen.add(path)
    text = path.read_text()
    result = {}
    for match in re.finditer(r'\\newlabel\{([^}]+)\}\{', text):
        body, _ = group(text, match.end() - 1)
        fields, pos = [], 0
        while pos < len(body) and body[pos] == '{':
            field, pos = group(body, pos)
            fields.append(field)
        result[match[1]] = fields
    for child in re.findall(r'\\@input\{([^}]+)\}', text):
        result.update(labels(root, root / child, seen))
    return result


def caption(root, item, book_labels):
    source = (root / item['source']).read_text()
    bodies = re.findall(r'\\begin\{figure\}.*?\\end\{figure\}', source, re.S)
    body = next(b for b in bodies if '\\label{' + item['id'] + '}' in b)
    start = body.index('\\caption{') + len('\\caption')
    text, _ = group(body, start)
    text = re.sub(r'\\(eqref|ref)\{([^}]+)\}',
                  lambda m: '(' + book_labels[m[2]][0] + ')' if m[1] == 'eqref'
                  else book_labels[m[2]][0], text)
    for macro, expansion in {'ZZ': r'\mathbb{Z}', 'FF': r'\mathbb{F}', 'Hh': r'\mathcal{H}'}.items():
        text = re.sub(r'\\' + macro + r'\b', lambda _: expansion, text)
    return ' '.join(text.replace('~', ' ').replace('--', '–').split())


def caption_html(text):
    pieces = re.split(r'(\$[^$]+\$)', text)
    return ''.join('<span class="figure-math">\\(' + html.escape(p[1:-1]) + '\\)</span>'
                   if p.startswith('$') else html.escape(p) for p in pieces)


def bundle(root, items, destination):
    files = {root / item['code'] for item in items}
    files.update(p.with_name('figure.tex') for p in list(files) if p.with_name('figure.tex').exists())
    files.update((root / 'listings').glob('*.py'))
    files.update(root / 'tools' / name for name in [
        'figure_style.py', 'applications_data.py', 'duality_data.py',
        'related_transform_data.py', 'tikz_figure.py'])
    files.add(root / 'data/flower.png')
    extra = {
        'requirements.txt': 'numpy\nscipy\nmatplotlib\nsympy\nPillow\npypdfium2\n',
        'README.md': '''# DAFT figure sources

The Discrete Algebra of the Fourier Transform, Gabriel Peyré.
An adapted English translation of L'algèbre discrète de la transformée de Fourier
(Ellipses, 2004).

Install Python 3.10 or later and run `python -m pip install -r requirements.txt`.
Run any `figures-code/chapter-*/figure-name/generate.py` with Python. Vector PDFs
are written to `figures/`; PNG previews are written to `build/figure-previews/`.
Each generator contains the parameters of the corresponding book figure.

The eight TikZ diagrams additionally require pdfLaTeX with the standalone,
New PX, and TikZ packages. Each diagram includes its `figure.tex` source.

The flower photograph in `data/flower.png` is the supplied input used for the
book's image-processing examples. Shared numerical routines are in `listings/`;
plotting conventions and common inputs are in `tools/`.

Book and gallery: https://mathematical-tours.github.io/daft/
''',
    }
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        contents = {p.relative_to(root).as_posix(): p.read_bytes() for p in files}
        contents.update({name: value.encode() for name, value in extra.items()})
        for name, data in sorted(contents.items()):
            info = zipfile.ZipInfo('daft-figures/' + name, (2026, 9, 9, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data)
    data = stream.getvalue()
    digest = hashlib.sha256(data).hexdigest()
    path = 'daft-sources/downloads/daft-figures-' + digest[:12] + '.zip'
    write(destination / path, data)
    return path, digest


def notebook(item, source, archive_path, archive_sha):
    setup = f'''from pathlib import Path
import hashlib
import importlib.util
import io
import os
import subprocess
import sys
import urllib.request
import zipfile

packages = {{"numpy": "numpy", "scipy": "scipy", "matplotlib": "matplotlib",
            "PIL": "Pillow", "sympy": "sympy"}}
missing = [package for module, package in packages.items()
           if importlib.util.find_spec(module) is None]
if missing:
    subprocess.check_call([sys.executable, "-m", "pip", "install", *missing])

archive_url = os.environ.get("DAFT_FIGURE_ARCHIVE", {str(PUBLIC + '/' + archive_path)!r})
archive_sha = {archive_sha!r}
workspace = Path(os.environ.get("DAFT_FIGURE_WORKDIR", "daft-figure-files")) / archive_sha[:12]
cache = workspace / "sources.zip"
workspace.mkdir(parents=True, exist_ok=True)
data = cache.read_bytes() if cache.exists() else urllib.request.urlopen(archive_url, timeout=60).read()
if hashlib.sha256(data).hexdigest() != archive_sha:
    raise ValueError("The figure-source download is incomplete or has changed. Download it again.")
cache.write_bytes(data)
with zipfile.ZipFile(io.BytesIO(data)) as archive:
    archive.extractall(workspace)
ROOT = (workspace / "daft-figures").resolve()
__file__ = str(ROOT / {item['code']!r})
print("Figure sources ready.")
'''
    output = f'''from IPython.display import Image, FileLink, display

display(Image(filename=str(ROOT / {item['preview_file']!r})))
display(FileLink(os.path.relpath(ROOT / {item['output']!r}, Path.cwd()),
                 result_html_prefix="Download the vector PDF: "))
'''
    cells = [
        ('markdown', f"# Figure {item['number']}: {item['title']}\n\n"
         f"*The Discrete Algebra of the Fourier Transform*, Gabriel Peyré.\n\n{item['caption']}\n\n"
         f"[Read this figure in the book]({PUBLIC}/daft-sources/DAFT-EN.pdf#page={item['pdf_page']})."),
        ('markdown', '## Prepare the example\n\nRun the cells in order in Jupyter or Google Colab. '
         'The first cell installs missing Python packages and downloads the figure sources and input data once. '
         'Subsequent runs use the cached download. Python 3.10 or later is required.'),
        ('code', setup),
        ('markdown', '## Generate the figure\n\nThis is the plotting code used by the book. '
         'Change its parameters and rerun this cell, then the display cell below, to explore the example.'),
        ('code', source),
        ('markdown', '## Display and save'),
        ('code', output),
    ]
    return {'nbformat': 4, 'nbformat_minor': 5,
            'metadata': {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                         'language_info': {'name': 'python', 'version': '3.12.0'},
                         'daft': {'figure': item['number'], 'source_sha256': hashlib.sha256(source.encode()).hexdigest()}},
            'cells': [dict(cell_type=kind, id=f'cell-{i}', metadata={}, source=text.splitlines(keepends=True),
                           **({'execution_count': None, 'outputs': []} if kind == 'code' else {}))
                      for i, (kind, text) in enumerate(cells)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--book-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--site-root', type=Path)
    args = parser.parse_args()
    root = args.book_root.resolve()
    dest = args.site_root.resolve() if args.site_root else root.parent
    items = json.loads((root / 'docs/figure-manifest.json').read_text())
    book_labels = labels(root, root / 'DAFT-EN.aux')
    book = PdfReader(root / 'DAFT-EN.pdf')
    archive_path, archive_sha = bundle(root, items, dest)
    figures = []
    for item in items:
        number, printed_page, _, anchor = book_labels[item['id']][:4]
        pdf_page = book.get_destination_page_number(book.named_destinations[anchor]) + 1
        if book.page_labels[pdf_page - 1] != printed_page:
            raise ValueError('Book pagination disagrees with the figure label: ' + item['id'])
        chapter, name = Path(item['output']).parts[1:3]
        title = TITLES[name]
        source = root / item['code']
        tikz = source.with_name('figure.tex')
        text = caption(root, item, book_labels)
        preview_path = f'img/daft-figures/{chapter}/{name}.webp'
        thumbnail_path = f'img/daft-figures/{chapter}/{name}-thumb.webp'
        with pdfium.PdfDocument(str(root / item['output'])) as doc:
            if len(doc) != 1:
                raise ValueError('A figure must have exactly one PDF page: ' + item['output'])
            page = doc[0]
            scale = min(2000 / page.get_width(), 1600 / page.get_height())
            bitmap = page.render(scale=scale)
            image = bitmap.to_pil().convert('RGB')
            (dest / preview_path).parent.mkdir(parents=True, exist_ok=True)
            image.save(dest / preview_path, 'WEBP', quality=92, method=6)
            width, height = image.size
            image.thumbnail((760, 560))
            image.save(dest / thumbnail_path, 'WEBP', quality=86, method=6)
            bitmap.close()
            page.close()
        item.update(number=number, pdf_page=pdf_page, title=title, caption=text,
                    preview_file=f'build/figure-previews/{chapter}/{name}.png')
        notebook_path = None
        if not tikz.exists():
            notebook_path = f'daft-sources/notebooks/{chapter}/{name}.ipynb'
            write(dest / notebook_path, json.dumps(notebook(item, source.read_text(), archive_path, archive_sha),
                                                   ensure_ascii=False, indent=1) + '\n')
        figures.append({
            'id': item['id'], 'number': number, 'chapter': int(number.split('.')[0]),
            'chapter_title': CHAPTERS[int(number.split('.')[0])], 'title': title,
            'caption': text, 'caption_html': caption_html(text),
            'book_page': printed_page, 'pdf_page': pdf_page,
            'pdf': '/daft-sources/' + item['output'],
            'preview': '/' + preview_path, 'thumbnail': '/' + thumbnail_path,
            'width': width, 'height': height, 'type': 'TikZ' if tikz.exists() else 'Python',
            'code': GITHUB + item['code'],
            'tikz': GITHUB + tikz.relative_to(root).as_posix() if tikz.exists() else None,
            'notebook': '/' + notebook_path if notebook_path else None,
            'colab': 'https://colab.research.google.com/github/mathematical-tours/mathematical-tours.github.io/blob/master/' + notebook_path if notebook_path else None,
            'pdf_sha256': hashlib.sha256((root / item['output']).read_bytes()).hexdigest(),
        })
    figures.sort(key=lambda f: tuple(map(int, f['number'].split('.'))))
    if len({f['number'] for f in figures}) != len(figures):
        raise ValueError('Duplicate figure numbers')
    data = {'book_sha256': hashlib.sha256((root / 'DAFT-EN.pdf').read_bytes()).hexdigest(),
            'archive': '/' + archive_path, 'archive_sha256': archive_sha,
            'chapters': [{'number': n, 'title': title, 'count': sum(f['chapter'] == n for f in figures)}
                         for n, title in CHAPTERS.items()], 'figures': figures}
    write(dest / '_data/daft_figures.json', json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    print(f'{len(figures)} figures; {sum(f["notebook"] is not None for f in figures)} notebooks; '
          f'{sum(f["type"] == "TikZ" for f in figures)} TikZ diagrams.')


if __name__ == '__main__':
    main()
