#!/usr/bin/env python3
"""One-time reversible migration of reading-edition references to preserved assets."""
from pathlib import Path
import hashlib,json,shutil
from collections import defaultdict
ROOT=Path(__file__).resolve().parents[1]

def main():
 project=ROOT/'figure-processing/figure-project.json';items=json.loads(project.read_text());stamp=ROOT/'figure-processing/migration-record.json'
 if stamp.exists():raise RuntimeError('Migration is already complete; use the regular book build.')
 by_source=defaultdict(list);changes=[]
 for i in items:
  for a in i['published_assets']:
   p=ROOT/'figures'/a['published'];assert hashlib.sha256(p.read_bytes()).hexdigest()==a['sha256'],p
  by_source[i['source']].append(i)
 # Compute every edit before changing any source.
 updates={}
 for source,figs in by_source.items():
  text=(ROOT/source).read_text();original=text
  for i in sorted(figs,key=lambda f:f['start'],reverse=True):
   if i['kind']=='cover':
    old=r'\includegraphics[width=\linewidth]{wave}';new=r'\includegraphics[width=\linewidth]{'+i['published_assets'][0]['published']+'}';assert text.count(old)==1;text=text.replace(old,new);continue
   segment=text[i['start']:i['end']];previous=segment
   for a,mapping in zip(i['assets'],i['published_assets']):
    if a['command']=='image':
     old=r'\image{'+a['args'][0]+'}{'+a['args'][1]+'}{'+a['args'][2]+'}'
     new=r'\includegraphics[width='+a['args'][1]+r'\linewidth]{'+mapping['published']+'}'
    else:old='{'+a['reference']+'}';new='{'+mapping['published']+'}'
    assert old in segment,(source,i['id'],old);segment=segment.replace(old,new)
   text=text[:i['start']]+segment+text[i['end']:]
   changes.append({'id':i['id'],'source':source,'asset_references':len(i['assets'])})
  updates[source]=text
 archive=ROOT/'figure-processing/archive';archive.mkdir(exist_ok=True)
 for name in ['figure-comparisons.pdf','figure-index.json','figure-comparisons.tex','comparison-layout.tex','README.md']:
  src=ROOT/'figure-processing'/name
  if src.exists():shutil.copy2(src,archive/('previous-'+name))
 legacy=archive/'legacy-assets';legacy.mkdir(exist_ok=True);moved=[]
 for p in sorted((ROOT/'figures').iterdir()):
  if p.name.startswith('chapter-'):continue
  if p.name=='.DS_Store':continue
  destination=legacy/p.name;assert not destination.exists();p.rename(destination);moved.append(str(destination.relative_to(ROOT)))
 for source,text in updates.items():(ROOT/source).write_text(text)
 p=ROOT/'book-preamble.tex';s=p.read_text().replace(r'\graphicspath{{./figures/}{./figures/ot/}}',r'\graphicspath{{./figures/}{./figure-processing/archive/legacy-assets/}{./figure-processing/archive/legacy-assets/ot/}}');p.write_text(s)
 # Correct a few historical folder misspellings without changing figure IDs.
 rename={'haard-scaling-wav':'haar-scaling-wavelet','ourier-approx-1d':'fourier-approximation-1d','iltering-optimal-curve':'filtering-optimal-curves','embeded-spaces-1d':'nested-spaces-1d','embeded-spaces-2d':'nested-spaces-2d'}
 substitutions={}
 for i in items:
  slug=i['figure_slug']
  if slug in rename:
   old=i['folder'];new=str(Path(old).with_name(rename[slug]));substitutions[old]=new
 for old,new in substitutions.items():
  for base in ['figures','figures-code']:(ROOT/base/old).rename(ROOT/base/new)
 # Paths can occur in reading sources, composite drivers and metadata.
 paths=[ROOT/'FundationsDataScience.tex',*[ROOT/s for s in by_source],*list((ROOT/'figures-code').rglob('*.tex'))]
 for p in dict.fromkeys(paths):
  text=p.read_text()
  for old,new in substitutions.items():text=text.replace(old+'/',new+'/')
  p.write_text(text)
 def update(value):
  if isinstance(value,str):
   for old,new in substitutions.items():value=value.replace(old,new)
   return value
  if isinstance(value,list):return [update(v) for v in value]
  if isinstance(value,dict):return {k:update(v) for k,v in value.items()}
  return value
 items=update(items)
 for i in items:
  if i['figure_slug'] in rename:i['figure_slug']=rename[i['figure_slug']]
  if not i['book_label']:i['caption_latex']=i['caption_source']
  i['in_book']=False
 project.write_text(json.dumps(items,indent=2,ensure_ascii=False)+'\n')
 stamp.write_text(json.dumps({'changes':changes,'archived_legacy_paths':moved,'renamed_folders':substitutions,'reading_edition':'preserved original assets; proposals excluded'},indent=2)+'\n')
 print('Migrated',sum(c['asset_references'] for c in changes)+1,'reading asset references;',len(moved),'legacy entries archived.')
if __name__=='__main__':main()
