#!/usr/bin/env python3
"""Audit figure coverage, preserved reading inputs and numerical/PDF artifacts."""
from pathlib import Path
import json,hashlib,re
import numpy as np
from pypdf import PdfReader
import pypdfium2 as pdf
ROOT=Path(__file__).resolve().parents[1];BUILD=ROOT/'build/figure-regeneration'

def unembedded_fonts(reader):
 """Inspect fonts in page resources and nested figure forms."""
 seen=set();missing=set()
 def resources(res):
  if res is None:return
  res=res.get_object()
  for ref in res.get('/Font',{}).values():
   font=ref.get_object()
   for f in font.get('/DescendantFonts',[font]):
    f=f.get_object()
    if f.get('/Subtype')=='/Type3':continue
    descriptor=f.get('/FontDescriptor')
    if not descriptor or not any(k in descriptor.get_object() for k in ('/FontFile','/FontFile2','/FontFile3')):
     missing.add(str(f.get('/BaseFont')))
  for ref in res.get('/XObject',{}).values():
   key=(getattr(ref,'idnum',None),getattr(ref,'generation',None))
   if key in seen:continue
   seen.add(key);resources(ref.get_object().get('/Resources'))
 for page in reader.pages:resources(page.get('/Resources'))
 return sorted(missing)

def publication_checks(issues):
 from build_book import chapter_list
 from figure_tex import braced
 book=PdfReader(ROOT/'FundationsDataScience.pdf');destinations=book.named_destinations
 compact=lambda text:re.sub(r'\s+','',text)
 date_pattern=r'(?:January|February|March|April|May|June|July|August|September|October|November|December)\d{1,2},\d{4}'
 date=re.search(date_pattern,compact(book.pages[0].extract_text())).group()
 paths=[ROOT/'FundationsDataScience.pdf',*[ROOT/'chapters-pdf'/f'{name}.pdf' for name in chapter_list()],ROOT/'docs/figures/figure-comparisons.pdf']
 results=[]
 for path in paths:
  reader=book if path.name=='FundationsDataScience.pdf' else PdfReader(path)
  names=reader.named_destinations;internal=external=0;dated=[]
  fonts=unembedded_fonts(reader)
  if fonts:issues.append(f'Unembedded fonts in {path.name}: {fonts}')
  for n,page in enumerate(reader.pages,1):
   content=page.extract_text()
   if re.search(r'\?\?|\[\?\]|Formula is false|\bTODO\b|\bToDo\b',content):issues.append(f'Unresolved text in {path.name}, page {n}')
   if date in compact(content):dated.append(n)
   for ref in page.get('/Annots',[]):
    annotation=ref.get_object();action=annotation.get('/A',{})
    if hasattr(action,'get_object'):action=action.get_object()
    target=action.get('/D',annotation.get('/Dest'))
    if action.get('/S')=='/GoToR':
     external+=1
     if target not in destinations:issues.append(f'Broken book link in {path.name}: {target}')
     link_file=action.get('/F')
     if isinstance(link_file,dict):link_file=link_file.get('/F')
     if (path.parent/str(link_file)).resolve()!=(ROOT/'FundationsDataScience.pdf').resolve():
      issues.append(f'Wrong book file in {path.name}: {link_file}')
    elif action.get('/S')=='/GoTo' or '/Dest' in annotation:
     internal+=1
     if isinstance(target,str) and target not in names:issues.append(f'Broken internal link in {path.name}: {target}')
    elif action.get('/S')=='/URI' and 'FundationsDataScience.pdf' in str(action.get('/URI','')):
     issues.append(f'Malformed reading-edition URI in {path.name}: {action.get("/URI")}')
  if dated!=[1]:issues.append(f'Unexpected edition-date pages in {path.name}: {dated}')
  results.append({'file':str(path.relative_to(ROOT)),'pages':len(reader.pages),'internal_links':internal,'book_links':external,'unembedded_fonts':fonts,'date_pages':dated})
 def labels(path):
  result={}
  for line in path.read_text().splitlines():
   if line.startswith(r'\newlabel{'):
    key,pos=braced(line,len(r'\newlabel'));body,_=braced(line,pos);number,_=braced(body,0);result[key]=number
  return result
 count=0
 for name in chapter_list():
  full=labels(ROOT/'build/book/chapters-tex'/f'{name}.aux');single=labels(ROOT/'build/chapters'/name/f'{name}.aux')
  for key,number in full.items():
   count+=1
   if single.get(key)!=number:issues.append(f'Chapter label mismatch {name}: {key}')
 builds=json.loads((ROOT/'build/build-report.json').read_text())
 for build in builds:
  if build['diagnostics']:issues.append('Compilation diagnostics in '+build['name'])
 return {'documents':results,'shared_labels_checked':count,'build_diagnostics':sum(len(b['diagnostics']) for b in builds)}

def main():
 items=json.loads((ROOT/'docs/figures/figure-project.json').read_text());by_id={i['id']:i for i in items};issues=[];report={};originals=set();props=set();expected=set();accepted=set();pending=set()
 allowed={'data/flower.png','data/figure-inputs/mandrill.png','data/figure-inputs/felix-1919.jpg'}
 for i in items:
  folder=i['folder'];assert i['code_directory']=='figures-code/'+folder;assert i['asset_directory']=='figures/'+folder
  state=i.get('review_state','revision-requested')
  if i.get('in_book')!=(state=='accepted'):issues.append('Acceptance state mismatch '+i['id'])
  for a in i['published_assets']:
   p=ROOT/'figures'/a['published'];compiled=p.with_name(p.stem+'-eps-converted-to.pdf') if p.suffix=='.eps' else p;originals.add(compiled.resolve())
   if hashlib.sha256(p.read_bytes()).hexdigest()!=a['sha256']:issues.append('Changed historical asset '+str(p))
   # Merged displays retain their old panels in the pending combined figure.
   merged_accepted=state=='merged' and by_id[i['merged_into']]['review_state']=='accepted'
   if state not in ('accepted','removed') and not merged_accepted:expected.add(compiled.resolve())
  proposal=ROOT/i['proposed_pdf'];props.add(proposal.resolve())
  if state=='accepted':expected.add(proposal.resolve());accepted.add(proposal.resolve())
  if state=='revision-requested':pending.add(proposal.resolve())
  if state in ('keep-original','removed','merged'):continue
  if not proposal.exists():issues.append('Missing reconstruction '+i['id']);continue
  d=PdfReader(proposal)
  if len(d.pages)!=1:issues.append('Multipage reconstruction '+i['id'])
  fonts=unembedded_fonts(d)
  if fonts:issues.append('Unembedded reconstruction fonts '+i['id']+': '+str(fonts))
  if i['engine']=='python':
   meta=json.loads((ROOT/i['asset_directory']/'provenance.json').read_text())
   for src in meta['data_sources']:
    if hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()!=src['sha256']:issues.append('Stale data '+i['id'])
    if Path(src['path']).suffix.lower() in {'.png','.jpg','.jpeg'} and src['path'] not in allowed:issues.append('Unapproved image input '+i['id']+': '+src['path'])
   if not (ROOT/i['code_directory']/'generate.py').exists():issues.append('Missing Python source '+i['id'])
  else:
   if not (ROOT/i['proposed_source']).exists():issues.append('Missing TeX composite '+i['id'])
   for t in i['tikz_sources']:
    if not (ROOT/t['source']).exists():issues.append('Missing TeX source '+i['id'])
 # The compiled reading edition must exactly implement the author's decisions.
 fls=ROOT/'build/book/FundationsDataScience.fls';inputs=set()
 for line in fls.read_text().splitlines():
  if line.startswith('INPUT '):
   p=Path(line[6:]);p=(ROOT/p).resolve() if not p.is_absolute() else p.resolve()
   if p.suffix.lower() in ['.pdf','.png','.jpg','.jpeg'] and ROOT in p.parents:inputs.add(p)
 for p in expected-inputs:issues.append('Accepted reading asset not compiled '+str(p))
 for p in inputs-expected:issues.append('Unexpected reading image '+str(p))
 for p in pending&inputs:issues.append('Unaccepted revision compiled '+str(p))
 if not accepted.issubset(inputs):issues.append('Accepted reconstructions missing from reading edition')
 active=[i for i in items if i['review_state'] not in ('removed','merged')]
 report.update(inventory_displays=len(items),active_displays=len(active),numbered_figures=len(set(i['figure_number'] for i in active if i['figure_number'])),figure_folders=len(set(i['folder'] for i in items)),accepted_reconstructions=sum(i['review_state']=='accepted' for i in items),retained_tikz_displays=sum(i['review_state']=='retained-tikz' for i in items),retained_original_displays=sum(i['review_state']=='keep-original' for i in items),removed_displays=sum(i['review_state']=='removed' for i in items),merged_displays=sum(i['review_state']=='merged' for i in items),historical_asset_files=len(originals),compiled_reading_image_files=len(inputs),unaccepted_proposals_in_book=len(pending&inputs),allowed_reconstruction_images=sorted(allowed))
 current=ROOT/'build/book/FundationsDataScience.pdf'
 # Publication comparison indexing and link destinations.
 index_path=ROOT/'docs/figures/figure-index.json'
 if index_path.exists():
  idx=json.loads(index_path.read_text())
  review_items=[i for i in items if i.get('comparison_visible',i['engine']!='tikz')]
  if len(idx)==len(review_items):
   comp=PdfReader(ROOT/'docs/figures/figure-comparisons.pdf');book=PdfReader(current);dest=book.named_destinations;badlinks=[];links=0
   for n,page in enumerate(comp.pages):
    for ref in page.get('/Annots',[]):
     annotation=ref.get_object();action=annotation.get('/A',{})
     if action.get('/S')=='/GoToR':
      links+=1;d=action.get('/D');target=action.get('/F')
      if d not in dest:badlinks.append({'page':n+1,'destination':str(d)})
      if str(target)!='../../FundationsDataScience.pdf':issues.append('Wrong comparison book path '+str(target))
   if badlinks:issues.append('Invalid comparison destinations '+str(badlinks[:6]))
   if len(comp.pages)!=len(review_items)+1:issues.append('Unexpected comparison page count')
   for i,j in zip(review_items,idx):
    if i['id']!=j['id'] or i['figure_number']!=j['figure_number'] or i['caption_latex']!=j['caption_latex']:issues.append('Comparison metadata mismatch '+i['id'])
    if i['book_anchor']:
     page=comp.pages[j['comparison_pdf_page']-1]
     actions=[ref.get_object().get('/A',{}) for ref in page.get('/Annots',[])]
     if not any(a.get('/S')=='/GoToR' and a.get('/D')==i['book_anchor'] for a in actions):issues.append('Missing figure-heading book link '+i['id'])
   report.update(comparison_displays=len(review_items),excluded_completed_displays=len(items)-len(review_items),comparison_pages=len(comp.pages),comparison_book_links=links,invalid_comparison_links=len(badlinks))
  else:issues.append('Comparison index still belongs to previous audit')
 report['publication']=publication_checks(issues)
 report['issues']=issues;(BUILD/'project-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));return bool(issues)
if __name__=='__main__':raise SystemExit(main())
