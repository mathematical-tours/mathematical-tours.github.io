#!/usr/bin/env python3
"""Prepare the Jekyll page, layout and cover using the current PDF downloads.

The destination is explicit: preview in build/site, or use the parent site root.
Only book-basics.md, _layouts/book-basics.html and its cover are written.
"""
from pathlib import Path
import argparse
import html
import json
import re
import shutil
from pypdf import PdfReader
import pypdfium2 as pdfium

ROOT = Path(__file__).resolve().parents[1]
CHAPTERS = [
    ('01-information', 'Claude Shannon and data compression',
     'Claude Shannon et la compression des données',
     'Binary codes, entropy, Huffman coding, and the limits of lossless compression.'),
    ('02-images', 'Image processing', 'Le traitement d’images',
     'Pixels, quantization, filtering, edges, colors, and elementary image transformations.'),
    ('03-sampling', 'Sparsity, inverse problems and compressed sensing',
     'Parcimonie, problèmes inverses et acquisition comprimée',
     'Sparse representations, regularization, and recovery from incomplete measurements.'),
    ('04-transport', 'Optimal transport', 'Transport optimal',
     'Assignments, transport plans, Wasserstein distances, and applications to distributions.'),
    ('05-neural-networks', 'Neural networks', 'Réseaux de neurones',
     'Classification, backpropagation, generative models, and connections with transport.'),
]
CHAPTER_FIGURES = [
    'basics-en-trees',
    'basics-en-canaux-coul',
    'basics-en-inpainting',
    'basics-en-image-eq',
    'basics-en-discriminative',
]

def pdf_stats(path):
    return len(PdfReader(path).pages), f'{path.stat().st_size/1_000_000:.1f}'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    dest = args.output.resolve()
    (dest/'_layouts').mkdir(parents=True, exist_ok=True)
    (dest/'img').mkdir(exist_ok=True)
    values = {}
    for lang, name in [('EN', 'IntroImaging'), ('FR', 'IntroImaging-FR')]:
        values[lang+'_PAGES'], values[lang+'_MB'] = pdf_stats(ROOT/(name+'.pdf'))
    metadata = (ROOT/'chapters-tex/shared/book-metadata.tex').read_text()
    values['DATE_FR'], values['DATE_EN'] = re.findall(r'\\newcommand\{\\BookEditionDate\}\{([^}]+)', metadata)
    rows = []
    gallery = json.loads((ROOT.parent/'_data/basics_figures.json').read_text())
    illustrations = {figure['id']: figure for figure in gallery['figures']}
    for n, (name, title, french, description) in enumerate(CHAPTERS, 1):
        links = []
        for lang, label in [('en', 'English PDF'), ('fr', 'PDF français')]:
            pages, _ = pdf_stats(ROOT/'chapters-pdf'/lang/(name+'.pdf'))
            url = "{{ '/book-basics-sources/chapters-pdf/"+lang+'/'+name+".pdf' | relative_url }}"
            accessible = f'Chapter {n}: {title} — English PDF' if lang=='en' else f'Chapitre {n} : {french} — PDF français'
            links.append(f'<a href="{url}" hreflang="{lang}" lang="{lang}" aria-label="{html.escape(accessible, quote=True)}">{label}<span class="pages">{pages} p.</span></a>')
        figure = illustrations[CHAPTER_FIGURES[n-1]]
        assert figure['chapter'] == n
        assert (ROOT.parent/figure['thumbnail'].lstrip('/')).is_file()
        image = "{{ '"+figure['thumbnail']+"' | relative_url }}"
        gallery_url = "{{ '/book-basics/figures/' | relative_url }}?chapter="+str(n)
        rows.append(
            '      <li>\n'
            f'        <a class="chapter-visual" href="{gallery_url}" aria-label="Browse the figures for chapter {n}: {html.escape(title, quote=True)}">'
            f'<img src="{image}" alt="{html.escape(figure["title"], quote=True)}" width="{figure["width"]}" height="{figure["height"]}" loading="lazy" decoding="async">'
            '<span>View figures <span aria-hidden="true">↗</span></span></a>\n'
            f'        <div class="chapter-copy"><span class="chapter-number">{n:02d}</span><div>'
            f'<h3 class="chapter-title">{html.escape(title)}<span class="chapter-title-fr" lang="fr">{html.escape(french)}</span></h3>'
            f'<p class="chapter-description">{html.escape(description)}</p></div></div>\n'
            '        <div class="chapter-actions">'+''.join(links)+'</div>\n'
            '      </li>'
        )
    values['CHAPTERS'] = '\n'.join(rows)
    source = (ROOT/'website/book-basics.md.in').read_text()
    for key, value in values.items():
        source = source.replace('@@'+key+'@@', str(value))
    assert not re.search(r'@@[A-Z_]+@@', source)
    (dest/'book-basics.md').write_text(source, encoding='utf8')
    shutil.copy2(ROOT/'website/book-basics.html', dest/'_layouts/book-basics.html')
    doc = pdfium.PdfDocument(ROOT/'IntroImaging.pdf')
    doc[0].render(scale=1.5).to_pil().convert('RGB').save(dest/'img/book-basics-cover.webp', quality=92)
    print(f'Prepared bilingual presentation page in {dest}')

if __name__ == '__main__':
    main()
