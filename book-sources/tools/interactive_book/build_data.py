#!/usr/bin/env python3
"""Compute discrete parameter grids, approved fallbacks, and standalone panels.

The web edition is a static artifact. Every slider stops at a computed parameter
value; images are never cross-faded or relabeled to simulate a calculation.
"""
from numerics import *
import argparse,shutil,os,time,traceback,concurrent.futures
import experiments_signal as S
import experiments_images as I
import experiments_learning as L
import pypdfium2 as pdfium
WEB=ROOT/'myst'
UPSTREAM_HASH=hashlib.sha256(b''.join(p.read_bytes() for p in sorted((ROOT/'figures-code').glob('*.py')))+b''.join(p.read_bytes() for p in sorted((ROOT/'data').rglob('*')) if p.is_file() and p.name!='.DS_Store')).hexdigest()
FAMILIES={**{k:(S,S.FUNCTIONS[k]) for k in S.FUNCTIONS},**{k:(I,I.FUNCTIONS[k]) for k in I.FUNCTIONS},**{k:(L,L.FUNCTIONS[k]) for k in L.FUNCTIONS}}


def static_figure(f):
    out=WEB/'assets/figures';out.mkdir(parents=True,exist_ok=True);paths=f['current_assets'];images_=[]
    for name in paths:
        path=ROOT/name
        # Prefer the already approved PNG when it is the same numerical asset.
        if path.name=='proposed.pdf' and path.with_suffix('.png').exists():path=path.with_suffix('.png')
        if path.suffix.lower()=='.pdf':
            pdf=pdfium.PdfDocument(path);page=pdf[0];im=page.render(scale=min(3,1500/page.get_width())).to_pil().convert('RGB');page.close();pdf.close()
        else:im=Image.open(path).convert('RGB')
        if im.width>1600:im.thumbnail((1600,2200),Image.Resampling.LANCZOS)
        images_.append(im)
    # The approved original.pdf is a composed figure for retained TikZ panels.
    # Its layout is already audited; avoid reconstructing its rows by guessing.
    composed=ROOT/f['asset_directory']/'original.pdf'
    if len(images_)>1 and composed.exists():
        pdf=pdfium.PdfDocument(composed);page=pdf[0];im=page.render(scale=min(3,1700/page.get_width())).to_pil().convert('RGB');page.close();pdf.close()
    else:im=images_[0]
    dst=out/(f['id']+'.webp');im.save(dst,lossless=True,method=4)
    return '../../assets/figures/'+dst.name


def build_case(f,force=False):
    cfg=f['interaction'];family=cfg['family'];view=cfg['view'];module,fn=FAMILIES[family];cs=module.controls(family,view)
    # Materialize only scientifically relevant controls. A parameter unrelated to
    # a selected view must never appear to control that view.
    if family=='wavelet2d' and view in ['wavalgo-2d-step-fwd','wavalgo-2d-step-bwd','wavalgo-2d','wavelets-2d']:cs=[c for c in cs if c['id']!='level']
    if family=='wavelet2d' and view=='wavelets-2d':cs=[c for c in cs if c['id']!='image']
    if family=='wavelet1d' and view in ['haard-scaling-wav','wavelets','inversion-transpose','haar-analysis-synthesis-matrix']:cs=[c for c in cs if c['id']!='signal']
    if family=='wavelet1d' and view in ['inversion-transpose','haar-analysis-synthesis-matrix']:
        cs[0]['values']=[1,2];cs[0]['default']=1
    if family=='models' and view in ['approximation-class','signal-model-sobolev']:cs=[c for c in cs if c['id']!='image']
    if (family=='noise' and view in ['noise-example-1d','sparse-signal']) or (family=='gaussian' and view=='filtering-denoised') or (family=='multiplicative' and view=='multnoise-hist-unstab'):cs=[c for c in cs if c['id']!='image']
    # For views that only vary a distribution, expose its mean directly.
    if family=='poisson' and view=='poisson-distrib':cs=[c for c in cs if c['id']!='image']
    if family=='coding' and view=='codage-differences':cs=[c for c in cs if c['id']!='skew']
    if family=='threshold' and view=='stein-optimal-curve':cs=[c for c in cs if c['id']!='method']
    if family=='threshold' and view=='wavthresh-2d-hard-vs-soft':
        cs[1]['values']=['Soft','Hard']
        cs[1]['default']=0
    if family=='compression' and view=='compression':cs[0]=slider('rate','Quantization step T',[.025,.05,.075,.1,.15,.25,.4])
    path=DATA/(f['id']+'.json')
    signature=hashlib.sha256((''.join(Path(v).read_text() for v in [__file__,module.__file__,ROOT/'tools/interactive_book/numerics.py'])+json.dumps(cs)+UPSTREAM_HASH).encode()).hexdigest()
    if path.exists() and not force:
        old=json.loads(path.read_text())
        if old.get('signature')==signature:return len(old['states'])
    states={};base={'image':'Flower','signal':'Piecewise smooth','level':1}
    for indices in itertools.product(*[range(len(c['values'])) for c in cs]):
        params={**base,**{c['id']:c['values'][index] for c,index in zip(cs,indices)}}
        result=fn(view,params);result['parameters']={c['id']:params[c['id']] for c in cs};states[','.join(map(str,indices))]=result
    payload={'id':f['id'],'family':family,'view':view,'controls':cs,'states':states,'seed':SEED,'signature':signature,'execution':'Computed Python parameter grid; no image interpolation.'}
    # Serialization rejects NaN/Infinity. A missing state is a build error.
    path.write_text(json.dumps(payload,ensure_ascii=False,separators=(',',':'),allow_nan=False)+'\n')
    return len(states)


def panel(f,static):
    import html
    folder=WEB/'live/panels';folder.mkdir(exist_ok=True);num=f['figure_number'];title=('Figure '+num if num else 'Illustration')+' · '+f['chapter_title'];cfg={'id':f['id'],'number':num,'title':title,'static':static,'data':'../data/'+f['id']+'.json' if f['interaction'] else None,'bookLabel':f['book_label']}
    encoded=json.dumps(cfg,ensure_ascii=False).replace('<','\\u003c')
    text=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{html.escape(title)}</title><link rel="stylesheet" href="../book-live.css"></head>
<body><main id="experiment" aria-label="{html.escape(title)}"><div class="toolbar"><div class="tabs" role="tablist" aria-label="Figure view"><button id="explore-tab" role="tab" aria-controls="explore" aria-selected="true">Explore</button><button id="static-tab" role="tab" aria-controls="static" aria-selected="false">Book figure</button></div><button id="reset" class="reset" aria-label="Reset all parameters">Reset</button></div>
<section id="explore" role="tabpanel" aria-labelledby="explore-tab" hidden><div id="controls" class="controls"></div><div id="status" class="status" role="status">Preparing experiment…</div><div id="panels" class="panels"></div><div id="metrics" class="metrics" aria-live="polite"></div><p id="note" class="note"></p><p class="grid-note">Each slider position shows a computed parameter value.</p></section>
<section id="static" role="tabpanel" aria-labelledby="static-tab"><img src="{static}" alt="{html.escape(title)} — approved book figure" class="static-image"><p class="static-note">The figure from the PDF edition.</p></section>
<noscript><p>JavaScript is needed for the parameter controls. The book figure remains available above.</p></noscript></main>
<script id="figure-config" type="application/json">{encoded}</script><script type="module" src="../book-live.mjs"></script></body></html>'''
    (folder/(f['id']+'.html')).write_text(text)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--family');parser.add_argument('--id');parser.add_argument('--force',action='store_true');parser.add_argument('--assets-only',action='store_true');args=parser.parse_args();DATA.mkdir(parents=True,exist_ok=True)
    figures=json.loads((WEB/'live/figures.json').read_text());report=[];start=time.time()
    hash_file=WEB/'assets/figure-hashes.json';static_hashes=json.loads(hash_file.read_text()) if hash_file.exists() else {}
    for f in figures:
        if args.id and f['id']!=args.id:continue
        if args.family and (not f['interaction'] or f['interaction']['family']!=args.family):continue
        static='../../assets/figures/'+f['id']+'.webp'
        dependencies=[ROOT/name for name in f['current_assets']]
        dependencies.extend(p.with_suffix('.png') for p in list(dependencies) if p.name=='proposed.pdf' and p.with_suffix('.png').exists())
        composed=ROOT/f['asset_directory']/'original.pdf'
        if len(f['current_assets'])>1 and composed.exists():dependencies.append(composed)
        digest=hashlib.sha256(b''.join(p.read_bytes() for p in dependencies)).hexdigest()
        if not (WEB/'assets/figures'/(f['id']+'.webp')).exists() or static_hashes.get(f['id'])!=digest or args.assets_only:static=static_figure(f)
        static_hashes[f['id']]=digest
        count=build_case(f,args.force) if f['interaction'] and not args.assets_only else 0
        panel(f,static);report.append({'id':f['id'],'states':count,'interactive':bool(f['interaction'])});print(f'{f["id"]}: {count} states',flush=True)
    hash_file.parent.mkdir(parents=True,exist_ok=True);hash_file.write_text(json.dumps(static_hashes,indent=2)+'\n')
    if not args.family and not args.id:
        pdfs=WEB/'assets/pdfs';pdfs.mkdir(parents=True,exist_ok=True)
        (pdfs/'chapters-pdf').mkdir(exist_ok=True)
        shutil.copyfile(ROOT/'FundationsDataScience.pdf',pdfs/'FundationsDataScience.pdf')
        for path in sorted((ROOT/'chapters-pdf').glob('*.pdf')):
            shutil.copyfile(path,pdfs/'chapters-pdf'/path.name)
            # Retire only our earlier generated flat copies. The chapter PDFs
            # contain links to ../FundationsDataScience.pdf, so keep that layout.
            flat=pdfs/path.name
            if flat.exists():flat.unlink()
        shutil.copyfile(ROOT/'figures/chapter-book-cover/paris-rooftops/artwork.jpg',WEB/'assets/cover.jpg')
        shutil.copyfile(ROOT/'data/bird.wav',WEB/'live/data/bird.wav')
        used_images=set()
        for data_file in DATA.glob('*.json'):
            d=json.loads(data_file.read_text())
            for state in d.get('states',{}).values():
                for panel_ in state['panels']:
                    for key in ['src','background']:
                        if panel_.get(key,'').startswith('images/'):used_images.add(panel_[key])
        for image_ in (DATA/'images').glob('*.webp'):
            if 'images/'+image_.name not in used_images:image_.unlink()
        (WEB/'data-report.json').write_text(json.dumps({'figures':report,'elapsed_seconds':time.time()-start},indent=2)+'\n')
    print(f'Finished in {time.time()-start:.1f} seconds',flush=True)
if __name__=='__main__':main()
