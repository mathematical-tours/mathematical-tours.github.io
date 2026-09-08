#!/usr/bin/env python3
"""Verify remote PDF destinations in the web edition's download tree."""
from pathlib import Path
from urllib.parse import unquote
import json
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[2]
HTML = ROOT / 'myst/_build/html'
PDFS = HTML / 'assets/pdfs'


def main():
    allowed = {p.resolve() for p in PDFS.rglob('*.pdf')}
    readers = {}
    errors = []
    checked = 0
    for chapter in sorted((PDFS / 'chapters-pdf').glob('*.pdf')):
        for page in PdfReader(chapter).pages:
            for reference in page.get('/Annots', []):
                action = reference.get_object().get('/A', {})
                if action.get('/S') != '/GoToR':
                    continue
                checked += 1
                name = action.get('/F', '')
                name = name.get('/F', '') if isinstance(name, dict) else str(name)
                target = (chapter.parent / unquote(name)).resolve()
                if target not in allowed or not target.is_file():
                    errors.append(f'{chapter.name}: unavailable PDF {name}')
                    continue
                if target not in readers:
                    readers[target] = PdfReader(target)
                destination = action.get('/D')
                if isinstance(destination, str):
                    valid = destination in readers[target].named_destinations
                else:
                    valid = bool(destination) and isinstance(destination[0], int) and 0 <= destination[0] < len(readers[target].pages)
                if not valid:
                    errors.append(f'{chapter.name}: missing destination {destination}')
    if not checked:
        errors.append('No chapter-to-book PDF links were found.')
    info = json.loads((HTML / 'build-info.json').read_text())
    report = {'content_sha256': info['content_sha256'],
              'remote_pdf_links_checked': checked, 'errors': errors}
    (ROOT / 'build/interactive-book/pdf-links.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
