#!/usr/bin/env python3
"""Create a small Jekyll input tree from the explicit publication inventory."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def regular_file(root, path):
    if root not in path.resolve().parents or path.is_symlink() or not path.is_file():
        raise ValueError(f'Not a regular file inside the source tree: {path}')
    if any(p.is_symlink() for p in path.parents if p != root and root in p.parents):
        raise ValueError(f'Symbolic link in source path: {path}')
    return path.relative_to(root).as_posix()


def collect(root, inventory):
    root = root.resolve()
    files, assets = set(), {}
    for kind in ('templates', 'assets'):
        for pattern in inventory[kind]:
            matches = sorted(root.glob(pattern))
            if not matches:
                raise ValueError(f'Publication input is missing: {pattern}')
            for path in matches:
                name = regular_file(root, path)
                files.add(name)
                if kind == 'assets':
                    assets[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_name = inventory['reader_manifest']
    manifest = root / manifest_name
    regular_file(root, manifest)
    reader = json.loads(manifest.read_text())
    for name, digest in reader['files'].items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError(f'Unsafe interactive-book path: {name}')
        path = manifest.parent / relative
        name = regular_file(root, path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'Interactive publication differs from its validated manifest: {name}')
        files.add(name)
        assets[name] = digest
    files.add(manifest_name)
    assets[manifest_name] = hashlib.sha256(manifest.read_bytes()).hexdigest()
    size = sum((root / name).stat().st_size for name in files)
    if size > inventory['maximum_bytes'] or len(files) > inventory['maximum_files']:
        raise ValueError(f'Publication input exceeds its budget: {size} bytes, {len(files)} files')
    return files, assets, size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    root = args.source.resolve()
    inventory = json.loads((Path(__file__).parent / 'public-site.json').read_text())
    files, assets, size = collect(root, inventory)
    if args.destination.exists():
        raise SystemExit('The destination must be a new directory; existing files are never removed.')
    args.destination.mkdir(parents=True)
    for name in sorted(files):
        target = args.destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / name, target)
    report = {'source_files': len(files), 'source_bytes': size, 'assets': assets,
              'maximum_bytes': inventory['maximum_bytes'], 'maximum_files': inventory['maximum_files']}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(f'Prepared {len(files)} publication inputs ({size / 1e6:.1f} MB).')


if __name__ == '__main__':
    main()
