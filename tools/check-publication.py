#!/usr/bin/env python3
"""Check the actual Pages payload, including every interactive asset and PDF."""
import argparse
import hashlib
import json
from pathlib import Path


def check(site, source_map):
    errors, files, size = [], [], 0
    for path in site.rglob('*'):
        if path.is_symlink():
            errors.append(f'Symbolic link in publication: {path.relative_to(site)}')
        elif path.is_file():
            files.append(path)
            size += path.stat().st_size
            if path.stat().st_nlink != 1:
                errors.append(f'Hard link in publication: {path.relative_to(site)}')
            if path.suffix in ('.tex', '.aux', '.log', '.fls', '.fdb_latexmk', '.toc', '.m', '.key', '.py'):
                errors.append(f'Development file in publication: {path.relative_to(site)}')
    for name, digest in source_map['assets'].items():
        path = site / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            errors.append(f'Missing or changed published asset: {name}')
    if size > source_map['maximum_bytes']:
        errors.append(f'Published size {size} exceeds the {source_map["maximum_bytes"]}-byte budget')
    if len(files) > source_map['maximum_files']:
        errors.append(f'Published file count {len(files)} exceeds the {source_map["maximum_files"]}-file budget')
    return {'files': len(files), 'bytes': size, 'verified_assets': len(source_map['assets']), 'errors': errors}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('site', type=Path)
    parser.add_argument('--source-map', type=Path, required=True)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    result = check(args.site.resolve(), json.loads(args.source_map.read_text()))
    text = json.dumps(result, indent=2) + '\n'
    print(text, end='')
    if args.report:
        args.report.write_text(text)
    if result['errors']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
