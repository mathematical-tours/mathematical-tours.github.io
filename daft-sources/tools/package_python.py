#!/usr/bin/env python3
"""Package the executable examples with their tests and small input data."""
from pathlib import Path
import json
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    target = ROOT / 'downloads/daft-python.zip'
    target.parent.mkdir(exist_ok=True)
    files = [ROOT/'requirements.txt', ROOT/'docs/program-verification.json']
    for directory in ['listings', 'tests']:
        files += sorted((ROOT/directory).glob('*.py'))
    files += [ROOT/'tools/verify_programs.py', ROOT/'data/flower.png']
    manifest = json.loads((ROOT/'docs/listing-manifest.json').read_text())
    manifest = [{key: item[key] for key in ('id', 'python_file', 'functions')}
                for item in manifest]
    readme = '''# DAFT — Python examples

Gabriel Peyré. Adapted English translation of L'algèbre discrète de la
transformée de Fourier (Ellipses, 2004, ISBN 9782729818678).

Use Python 3.12. From this folder:

    python -m pip install -r requirements.txt
    python -W error tools/verify_programs.py

For an example, add listings/ to Python's module search path:

    import sys
    sys.path.insert(0, "listings")
    from fft_algorithms import fft_dit
    print(fft_dit([1, 2, 0, -1]))

All 88 printed programs are represented by these modules. Shared functions
and imports are included. The tests use independent mathematical references.
Full book and figure sources: https://mathematical-tours.github.io/daft/
'''
    with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        generated = {
            'README.md': readme,
            'docs/listing-manifest.json': json.dumps(manifest, indent=2)+'\n',
            'data/README.md': 'The flower photograph was supplied with the companion MFDS sources.\nNumerical examples convert it to grayscale and resize it.\n',
        }
        for name, content in generated.items():
            info = zipfile.ZipInfo('daft-python/'+name, (2026, 9, 8, 0, 0, 0))
            archive.writestr(info, content, compress_type=zipfile.ZIP_DEFLATED)
        for path in files:
            info = zipfile.ZipInfo('daft-python/'+path.relative_to(ROOT).as_posix(),
                                   (2026, 9, 8, 0, 0, 0))
            archive.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED)
    print(f'{target.name}: {len(files)+len(generated)} files, {target.stat().st_size:,} bytes')


if __name__ == '__main__':
    main()
