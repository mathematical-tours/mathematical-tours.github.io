#!/usr/bin/env python3
"""Release checks for coverage, parameter grids, local links, and PDF fidelity."""
from pathlib import Path
import json,sys,itertools,hashlib,collections,re,urllib.parse,html
from html.parser import HTMLParser
ROOT=Path(__file__).resolve().parents[2];WEB=ROOT/'myst';errors=[]

def check(ok,message):
    if not ok:errors.append(message)

class Page(HTMLParser):
    def __init__(self,path):super().__init__();self.ids=set();self.links=[];self.frames=[];self.resources=[];self.math_errors=[];self.feed(path.read_text())
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if 'id' in a:self.ids.add(a['id'])
        if tag=='a' and 'href' in a:self.links.append(a['href'])
        if tag=='iframe' and 'src' in a:self.frames.append(a['src'])
        if tag in ['script','img'] and a.get('src'):self.resources.append(a['src'])
        if tag=='link' and a.get('rel')=='stylesheet' and a.get('href'):self.resources.append(a['href'])
        if 'katex-error' in a.get('class',''):self.math_errors.append(a.get('title','Unknown math error'))

def main():
    figures=json.loads((WEB/'live/figures.json').read_text());conv=json.loads((WEB/'conversion-report.json').read_text());mapping=json.loads((WEB/'interactions.json').read_text())
    check(len(conv['chapters'])==19,'The active book must contain all 19 chapters.')
    check(not conv['unknown_tex'],'Unconverted LaTeX remains.')
    check(not conv['warnings'],'Pandoc emitted conversion warnings.')
    check(len(figures)==len(conv['figures']),'Figure count changed during conversion.')
    check(len(set(conv['figures']))==len(conv['figures']),'A figure was inserted twice.')
    check(set(conv['figures'])=={f['id'] for f in figures},'The reading edition dropped a figure.')
    checks=[];state_count=0
    check(set(mapping)=={f['id'] for f in figures if f['review_state']=='accepted'},'The interaction registry does not cover exactly the accepted numerical figures.')
    for f in figures:
        check((WEB/'assets/figures'/(f['id']+'.webp')).exists(),'Missing static fallback '+f['id'])
        if not f['interaction']:continue
        path=WEB/'live/data'/(f['id']+'.json');built=WEB/'_build/html/live/data'/path.name
        check(built.exists() and built.read_bytes()==path.read_bytes(),'Stale built experiment data: '+f['id'])
        d=json.loads(path.read_text());cs=d['controls'];states=d['states'];expected={','.join(map(str,v)) for v in itertools.product(*[range(len(c['values'])) for c in cs])}
        check(set(states)==expected,'Incomplete Cartesian parameter grid '+f['id'])
        check(all(0<=c['default']<len(c['values']) for c in cs),'Invalid default control '+f['id'])
        for key,state in states.items():
            for p in state['panels']:
                if p['kind']=='image':check((WEB/'live/data'/p['src']).exists(),'Missing image '+p['src'])
                if p['kind']=='mesh':check((WEB/'live/data'/p['background']).exists(),'Missing mesh background')
                if p['kind']=='plot':
                    for s in p['series']:check(len(s['x'])==len(s['y']),'Inconsistent curve lengths '+f['id'])
        for j,c in enumerate(cs):
            groups=collections.defaultdict(set)
            for k,state in states.items():
                v=k.split(',');v[j]='*';groups[','.join(v)].add(json.dumps(state['panels'],sort_keys=True))
            check(max(map(len,groups.values()))>1,f'Control has no visual effect: {f["id"]}/{c["id"]}')
        state_count+=len(states);checks.append({'id':f['id'],'figure_number':f['figure_number'],'family':d['family'],'states':len(states),'controls':[c['label'] for c in cs]})
    root=WEB/'_build/html';build_info=json.loads((root/'build-info.json').read_text());base=build_info.get('base_url','');pages={p:Page(p) for p in [root/'index.html',*root.glob('*/index.html')] if p.parent.name not in ['live','assets']}
    broken=[];frame_count=0
    for path,page in pages.items():
        check(not page.math_errors,'Math rendering error '+str(path))
        frame_count+=len(page.frames)
        for url in page.resources:
            check(not urllib.parse.urlsplit(url).scheme.startswith('http'),'External reading dependency: '+url)
        for url in [*page.links,*page.frames,*page.resources]:
            u=urllib.parse.urlsplit(html.unescape(url))
            if u.scheme or u.netloc or not u.path and not u.fragment:continue
            local_path=u.path[len(base):] if base and u.path.startswith(base+'/') else u.path
            target=(root/local_path.lstrip('/')) if local_path.startswith('/') else path.parent/local_path
            target=target.resolve()
            if target.is_dir():target=target/'index.html'
            if not target.exists() and not target.suffix:target=target/'index.html'
            if not target.exists():broken.append((str(path.relative_to(root)),url));continue
            if u.fragment and target.suffix=='.html' and '/live/' not in str(target):
                other=pages.get(target) or Page(target)
                if urllib.parse.unquote(u.fragment) not in other.ids:broken.append((str(path.relative_to(root)),url))
    check(not broken,'Broken local links: '+str(broken[:20]))
    check(frame_count==len(figures),'Built HTML iframe count does not match the book.')
    pdf_hashes={}
    for p in [ROOT/'FundationsDataScience.pdf',*sorted((ROOT/'chapters-pdf').glob('*.pdf'))]:
        digest=hashlib.sha256(p.read_bytes()).hexdigest();copy=root/'assets/pdfs'/p.relative_to(ROOT)
        check(copy.exists() and hashlib.sha256(copy.read_bytes()).hexdigest()==digest,'Published PDF changed or is absent: '+p.name);pdf_hashes[p.name]=digest
    report={'chapters':len(conv['chapters']),'numbered_equations':len(conv['equations']),'theorem_blocks':sum(c['theorems'] for c in conv['chapters']),'figure_displays':len(figures),'interactive_figures':len(checks),'static_figures':len(figures)-len(checks),'computed_states':state_count,'families':len({c['family'] for c in checks}),'broken_links':broken,'errors':errors,'pdf_sha256':pdf_hashes,'interactions':checks}
    report['content_sha256']=build_info['content_sha256']
    out=ROOT/'build/interactive-book';out.mkdir(exist_ok=True);(out/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['interactions','pdf_sha256']},indent=2))
    if errors:sys.exit(1)
if __name__=='__main__':main()
