#!/usr/bin/env python3
"""Check that the public DAFT gallery matches the book and its figure sources."""
import hashlib
import json
from pathlib import Path
import zipfile


def check(root):
    data = json.loads((root / '_data/daft_figures.json').read_text())
    manifest = json.loads((root / 'daft-sources/docs/figure-manifest.json').read_text())
    sources = {f['id']: f for f in manifest}
    figures = data['figures']
    errors = []
    def require(condition, message):
        if not condition:
            errors.append(message)
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None

    require(data['book_sha256'] == digest(root / 'daft-sources/DAFT-EN.pdf'),
            'The book changed: rebuild the gallery to refresh PDF page links.')
    require(set(sources) == {f['id'] for f in figures} and len(figures) == len(sources),
            'The gallery must contain every numbered figure exactly once.')
    numbers = [tuple(map(int, f['number'].split('.'))) for f in figures]
    require(numbers == sorted(set(numbers)), 'Figure numbers must be unique and in book order.')
    archive_path = root / data['archive'].lstrip('/')
    require(digest(archive_path) == data['archive_sha256'], 'The source archive checksum changed.')
    notebooks = 0
    with zipfile.ZipFile(archive_path) as archive:
        for figure in figures:
            item = sources[figure['id']]
            for key in ['pdf', 'preview', 'thumbnail']:
                require((root / figure[key].lstrip('/')).is_file(), 'Missing ' + figure[key])
            require(figure['pdf_sha256'] == digest(root / figure['pdf'].lstrip('/')),
                    'A figure changed: rebuild its gallery preview: ' + figure['number'])
            require(figure['chapter'] == int(figure['number'].split('.')[0]),
                    'Incorrect chapter for Figure ' + figure['number'])
            source = (root / 'daft-sources' / item['code']).read_text()
            require(archive.read('daft-figures/' + item['code']).decode() == source,
                    'The archived generator changed: ' + item['code'])
            if figure['type'] == 'Python':
                notebook = json.loads((root / figure['notebook'].lstrip('/')).read_text())
                code = ''.join(notebook['cells'][4]['source'])
                require(code == source, 'Notebook plotting code differs from the book: ' + figure['number'])
                require(data['archive_sha256'] in ''.join(notebook['cells'][2]['source']),
                        'Notebook points to another source archive: ' + figure['number'])
                require(notebook['metadata']['daft']['figure'] == figure['number'],
                        'Incorrect notebook figure number: ' + figure['number'])
                require(all(not c.get('outputs') for c in notebook['cells']),
                        'Notebook contains stored execution output: ' + figure['number'])
                notebooks += 1
            else:
                require(figure['type'] == 'TikZ' and figure['tikz'] and not figure['notebook'],
                        'Missing TikZ source link: ' + figure['number'])
                path = Path(item['code']).with_name('figure.tex')
                require(archive.read('daft-figures/' + path.as_posix()) ==
                        (root / 'daft-sources' / path).read_bytes(),
                        'The archived TikZ source changed: ' + figure['number'])
        for path in ['data/flower.png', 'tools/figure_style.py', 'tools/applications_data.py',
                     'tools/duality_data.py', 'tools/related_transform_data.py', 'tools/tikz_figure.py']:
            require(archive.read('daft-figures/' + path) == (root / 'daft-sources' / path).read_bytes(),
                    'The archived figure dependency changed: ' + path)
        for path in (root / 'daft-sources/listings').glob('*.py'):
            require(archive.read('daft-figures/listings/' + path.name) == path.read_bytes(),
                    'The archived numerical module changed: ' + path.name)
    return {'figures': len(figures), 'notebooks': notebooks, 'errors': errors}


if __name__ == '__main__':
    result = check(Path(__file__).resolve().parents[1])
    print(json.dumps(result, indent=2))
    raise SystemExit(bool(result['errors']))
