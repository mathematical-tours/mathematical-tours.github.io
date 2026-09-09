#!/usr/bin/env python3
"""Regenerate registered DAFT figures; select labels with --figure."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure',action='append',default=[])
    args=parser.parse_args()
    manifest=json.loads((ROOT/'docs/figure-manifest.json').read_text())
    selected=[f for f in manifest if f['status']=='regenerated' and
              (not args.figure or f['id'] in args.figure)]
    missing=set(args.figure)-{f['id'] for f in selected}
    if missing:
        parser.error('Unregistered figures: '+', '.join(sorted(missing)))
    for f in selected:
        subprocess.run([sys.executable,'-W','error',str(ROOT/f['code'])],
                       cwd=ROOT,check=True)
        if not (ROOT/f['output']).is_file():
            raise FileNotFoundError(f['output'])
        print(f['id'],flush=True)
    print(f'Regenerated {len(selected)} figures.')

if __name__=='__main__':
    main()
