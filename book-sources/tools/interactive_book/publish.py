#!/usr/bin/env python3
"""Copy the validated interactive book into the existing GitHub Pages tree."""
from pathlib import Path
import hashlib
import json
import shutil
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'myst/_build/html'
DESTINATION = ROOT.parent / 'book/interactive'
URL = 'https://mathematical-tours.github.io/book/interactive/'


def main():
    info = json.loads((SOURCE / 'build-info.json').read_text())
    prefix = urlsplit(URL).path.rstrip('/')
    if info['base_url'] != prefix:
        raise SystemExit(f'Rebuild with BASE_URL={prefix} before publication.')
    if info.get('public_url') != URL:
        raise SystemExit('The export metadata must use the public book URL.')
    audit = json.loads((ROOT / 'docs/interactive-book/validation.json').read_text())
    for report in [audit, audit['browser'], audit['reading'], audit['pdf_links']]:
        if (report.get('errors') or report.get('issues') or report.get('external_requests')
                or report.get('content_sha256') != info['content_sha256']):
            raise SystemExit('Publication requires successful audits of this exact build.')
        if 'book_url' in report and urlsplit(report['book_url']).path.rstrip('/') != prefix:
            raise SystemExit('Run the browser checks at the publication prefix first.')
    if audit['browser']['passed'] != audit['figure_displays']:
        raise SystemExit('The browser audit is incomplete.')

    # Keep publication ownership explicit. Only files listed by the previous
    # manifest may be replaced or pruned; unrelated files are never removed.
    manifest_path = DESTINATION / 'publication.json'
    old = json.loads(manifest_path.read_text())['files'] if manifest_path.exists() else {}
    if DESTINATION.is_symlink():
        raise SystemExit('The publication directory must not be a symbolic link.')
    for name in old:
        if Path(name).is_absolute() or '..' in Path(name).parts:
            raise SystemExit('Unsafe path in the previous publication manifest.')
    if DESTINATION.exists():
        unexpected = [p for p in DESTINATION.rglob('*')
                      if p.is_symlink() or p.is_file() and p != manifest_path
                      and p.relative_to(DESTINATION).as_posix() not in old]
        if unexpected:
            raise SystemExit(f'Unmanaged files in the publication directory: {unexpected[:3]}')

    files = {}
    for source in sorted(SOURCE.rglob('*')):
        if source.is_symlink():
            raise SystemExit(f'Symbolic link in the static export: {source}')
        if not source.is_file():
            continue
        name = source.relative_to(SOURCE).as_posix()
        files[name] = hashlib.sha256(source.read_bytes()).hexdigest()
    DESTINATION.mkdir(parents=True, exist_ok=True)
    for name in files:
        target = DESTINATION / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(SOURCE / name, target)
    for name in old.keys() - files.keys():
        (DESTINATION / name).unlink(missing_ok=True)
    for directory in sorted(DESTINATION.rglob('*'), reverse=True):
        if directory.is_dir() and not any(directory.iterdir()):
            directory.rmdir()
    for name, digest in files.items():
        if hashlib.sha256((DESTINATION / name).read_bytes()).hexdigest() != digest:
            raise SystemExit(f'Publication copy differs from the validated build: {name}')
    manifest = {**info, 'url': URL, 'files': files}
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Prepared {len(files)} verified files at {DESTINATION}')
    print(f'Commit and push the publication and source changes to publish {URL}')


if __name__ == '__main__':
    main()
