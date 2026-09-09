#!/usr/bin/env python3
"""Build the bilingual books and standalone chapters, with isolated TeX outputs."""
from pathlib import Path
import argparse
import concurrent.futures
import json
import os
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
CHAPTERS = ('01-information', '02-images', '03-sampling', '04-transport', '05-neural-networks')

def build(target):
    source, destination = target
    out = ROOT / 'build/latex' / source.stem
    for lang in ('en', 'fr'):
        (out / 'chapters-tex' / lang).mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['PATH'] = os.pathsep.join((env.get('PATH', ''), '/Library/TeX/texbin', '/opt/homebrew/bin'))
    command = ['latexmk', '-pdf', '-recorder', '-interaction=nonstopmode', '-halt-on-error',
               f'-outdir={out.relative_to(ROOT)}', str(source.relative_to(ROOT))]
    with (out/'build.log').open('w') as stream:
        subprocess.run(command, cwd=ROOT, env=env, stdout=stream,
                       stderr=subprocess.STDOUT, check=True)
    log = (out/(source.stem+'.log')).read_text(errors='replace')
    problems = re.findall(r'^(?:.*Warning:.*|Overfull .*|Underfull .*|Missing character:.*|!.*)$', log, re.M)
    if problems:
        raise RuntimeError(f'{source.name}: unresolved TeX diagnostics:\n'+'\n'.join(problems))
    blg = out/(source.stem+'.blg')
    if blg.exists():
        text = blg.read_text(errors='replace')
        if re.search(r'Warning--|error message', text):
            raise RuntimeError(f'Bibliography diagnostics in {blg}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(out/(source.stem+'.pdf'), destination)
    print(f'Built {destination.relative_to(ROOT)}', flush=True)
    return {'source':str(source.relative_to(ROOT)), 'pdf':str(destination.relative_to(ROOT)),
            'log':str((out/(source.stem+'.log')).relative_to(ROOT)), 'warnings':0}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--books', action='store_true', help='build only the complete books')
    selection.add_argument('--chapters', action='store_true', help='build only the standalone chapters')
    parser.add_argument('--jobs', type=int, default=2)
    args=parser.parse_args()
    targets=[]
    if not args.chapters:
        targets.extend((ROOT/(s+'.tex'),ROOT/(s+'.pdf')) for s in ('IntroImaging','IntroImaging-FR'))
    if not args.books:
        targets.extend((ROOT/'chapters-tex/standalone'/f'{lang}-{name}.tex',
                        ROOT/'chapters-pdf'/lang/f'{name}.pdf')
                       for lang in ('en','fr') for name in CHAPTERS)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results=list(pool.map(build,targets))
    for old, current in json.loads((ROOT/'tools/legacy-downloads.json').read_text()).items():
        source = ROOT/current
        if source.exists():
            destination = ROOT/old
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    (ROOT/'reports').mkdir(exist_ok=True)
    report='build-books.json' if args.books else 'build-chapters.json' if args.chapters else 'build.json'
    (ROOT/'reports'/report).write_text(json.dumps(results,indent=2)+'\n')

if __name__=='__main__':
    main()
