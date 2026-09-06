#!/usr/bin/env python3
"""Regenerate accepted and pending figures from their per-figure Python/TeX sources."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT/'figure-processing/figure-project.json'
BUILD=ROOT/'build/figure-regeneration'
DIAGNOSTIC=re.compile(r'^(?:LaTeX(?: Font)? Warning:|Package .+ Warning:|pdfTeX warning|Overfull \\[hv]box|Underfull \\[hv]box|Missing character:)')

def compile_tex(source,output,job):
    directory=BUILD/'tex'/job;directory.mkdir(parents=True,exist_ok=True)
    source=ROOT/source;output=ROOT/output
    from figure_inventory import graphics
    dependencies=[source,ROOT/'book-preamble.tex',ROOT/'figures-code/tikz-preamble.tex']
    dependencies.extend(ROOT/a['path'] for a in graphics(source.read_text()) if a['path'])
    digest=hashlib.sha256(b''.join(p.read_bytes() for p in dependencies)).hexdigest()
    stamp=directory/(source.stem+'.sha256')
    if output.exists() and stamp.exists() and stamp.read_text()==digest:return []
    command=['pdflatex','-interaction=nonstopmode','-halt-on-error','-file-line-error','-recorder',
             f'-output-directory={directory}',str(source)]
    env=os.environ.copy();env['max_print_line']='1000'
    for i in range(2):
        p=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True)
        (directory/(source.stem+f'-{i}.txt')).write_text(p.stdout+p.stderr)
        if p.returncode:raise RuntimeError('\n'.join(p.stdout.splitlines()[-16:]))
    log=(directory/(source.stem+'.log')).read_text(errors='replace')
    diagnostics=[line for line in log.splitlines() if DIAGNOSTIC.match(line)]
    if diagnostics:raise RuntimeError(str(source)+'\n'+'\n'.join(diagnostics))
    shutil.copy2(directory/(source.stem+'.pdf'),output);stamp.write_text(digest)
    return diagnostics

def build_entry(item,originals_only=False):
    compile_tex(item['original_source'],item['original_pdf'],item['id'])
    if originals_only:return 'original'
    if item.get('review_state') in ('removed','merged','keep-original'):return 'retained-original'
    code=ROOT/item['code_directory']
    generator=code/'generate.py'
    if generator.exists():
        provenance=ROOT/item['asset_directory']/'provenance.json'
        dependencies=[generator,ROOT/'figures-code/requirements.txt',*sorted((ROOT/'figures-code').glob('*.py'))]
        if provenance.exists():
            meta=json.loads(provenance.read_text());dependencies.extend(ROOT/d['path'] for d in meta.get('data_sources',[]))
        # A changed generator may deliberately stop using a previous data file.
        # Missing cached dependencies invalidate the cache; the generator then
        # records its current inputs, or reports a real missing-input error.
        valid_dependencies=all(p.exists() for p in dependencies)
        digest=hashlib.sha256(b''.join(p.read_bytes() for p in dependencies if p.exists())).hexdigest()
        stamp=BUILD/(item['id']+'-python.sha256')
        outputs=[ROOT/item['proposed_pdf'],provenance,ROOT/item['asset_directory']/'proposed.png']
        if valid_dependencies and all(p.exists() for p in outputs) and stamp.exists() and stamp.read_text()==digest:return 'python'
        env=os.environ.copy();env.update(MPLCONFIGDIR=str(BUILD/'matplotlib'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        p=subprocess.run([sys.executable,str(generator)],cwd=ROOT,env=env,capture_output=True,text=True)
        (BUILD/(item['id']+'-python.txt')).write_text(p.stdout+p.stderr)
        if p.returncode or re.search(r'(?:Warning:|findfont:)',p.stderr):raise RuntimeError(p.stdout+p.stderr)
        if not (ROOT/item['proposed_pdf']).exists():raise RuntimeError('Generator did not create expected PDF')
        # Recompute with newly recorded data dependencies on the first generation.
        meta=json.loads(provenance.read_text())
        dependencies=[generator,ROOT/'figures-code/requirements.txt',*sorted((ROOT/'figures-code').glob('*.py')),*[ROOT/d['path'] for d in meta.get('data_sources',[])]]
        stamp.write_text(hashlib.sha256(b''.join(p.read_bytes() for p in dependencies)).hexdigest())
        return 'python'
    if item['engine']=='tikz':
        for entry in item['tikz_sources']:compile_tex(entry['source'],entry['output'],item['id'])
        compile_tex(item['proposed_source'],item['proposed_pdf'],item['id'])
        return 'tikz'
    return 'pending'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chapter');parser.add_argument('--id');parser.add_argument('--jobs',type=int,default=4)
    parser.add_argument('--originals-only',action='store_true');parser.add_argument('--allow-pending',action='store_true')
    args=parser.parse_args();items=json.loads(PROJECT.read_text())
    selected=[i for i in items if (not args.chapter or i['chapter'] in args.chapter.split(',')) and (not args.id or i['id']==args.id)]
    results=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures={pool.submit(build_entry,i,args.originals_only):i for i in selected}
        for future in as_completed(futures):
            i=futures[future]
            try:
                engine=future.result();result={'id':i['id'],'number':i['figure_number'],'engine':engine,'issues':[]}
                if engine in ('python','tikz'):
                    i['engine']=engine;i['status']=i.get('review_state','generated-awaiting-audit')
                    provenance=ROOT/i['asset_directory']/'provenance.json'
                    if provenance.exists():
                        meta=json.loads(provenance.read_text());i['notes']=meta.get('notes',i['notes']);i['data_sources']=meta.get('data_sources',[])
                print(i['figure_number'] or 'Unnumbered',i['id'],engine,flush=True)
            except Exception as error:
                result={'id':i['id'],'number':i['figure_number'],'engine':'failed','issues':[str(error)]}
                print('FAILED',i['id'],str(error),flush=True)
            results.append(result)
    fresh=json.loads(PROJECT.read_text());updated={i['id']:i for i in selected}
    for item in fresh:
        if item['id'] in updated:
            for field in ('engine','status','notes','data_sources'):
                if field in updated[item['id']]:item[field]=updated[item['id']][field]
    PROJECT.write_text(json.dumps(fresh,indent=2,ensure_ascii=False)+'\n')
    report=BUILD/('build-'+(args.id or args.chapter or 'all')+'.json');report.write_text(json.dumps(results,indent=2)+'\n')
    return 1 if any(r['issues'] or (r['engine']=='pending' and not args.allow_pending) for r in results) else 0

if __name__=='__main__':raise SystemExit(main())
