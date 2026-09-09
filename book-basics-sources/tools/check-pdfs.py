#!/usr/bin/env python3
"""Check PDF metadata, text, embedded fonts, figures and standalone numbering."""
from pathlib import Path
import json,re,unicodedata
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[1]
date_fr, date_en = re.findall(r'\\newcommand\{\\BookEditionDate\}\{([^}]+)',
                            (ROOT/'chapters-tex/shared/book-metadata.tex').read_text())

def embedded_fonts(resources,seen=None):
 seen=set() if seen is None else seen
 found=[]
 if not resources:return found
 resources=resources.get_object()
 for font in resources.get('/Font',{}).get_object().values() if '/Font' in resources else []:
  f=font.get_object();base=str(f.get('/BaseFont','unknown'));sub=str(f.get('/Subtype'))
  if sub=='/Type3': found.append((base,'Type3'));continue
  descendants=f.get('/DescendantFonts',[f])
  for d in descendants:
   desc=d.get_object().get('/FontDescriptor')
   ok=bool(desc and any(k in desc.get_object() for k in ('/FontFile','/FontFile2','/FontFile3')))
   found.append((base,'embedded' if ok else 'NOT EMBEDDED'))
 for obj in resources.get('/XObject',{}).get_object().values() if '/XObject' in resources else []:
  if hasattr(obj,'idnum'):
   ident=(obj.idnum,obj.generation)
   if ident in seen:continue
   seen.add(ident)
  x=obj.get_object()
  if x.get('/Subtype')=='/Form':found+=embedded_fonts(x.get('/Resources'),seen)
 return found

records=[]
for source in list(ROOT.glob('IntroImaging*.pdf'))+sorted((ROOT/'chapters-pdf').rglob('*.pdf')):
 lang='fr' if source.name=='IntroImaging-FR.pdf' or source.parent.name=='fr' else 'en'
 reader=PdfReader(source);texts=[p.extract_text() or '' for p in reader.pages]
 text='\n'.join(texts)
 assert all('�' not in s for s in texts),source
 assert not any(bad in text for bad in ('donnÃ','Ã©','Ã¨','??')),source
 assert 'Gabriel Peyré' in texts[0],source
 assert reader.metadata.author=='Gabriel Peyré',source
 date=date_fr if lang=='fr' else date_en
 assert text.count(date)==1,(source,'date count',text.count(date))
 if lang=='fr':assert 'sciences des données' in texts[0],source
 if source.parent.name in ('en','fr'):
  number=int(source.name[:2]);label='CHAPITRE' if lang=='fr' else 'CHAPTER'
  assert f'{label} {number}' in texts[0],(source,'chapter number')
  assert 'CNRS & DMA' in texts[0],source
 else:
  for n in range(1,6):assert f'CHAPITRE {n}' in text or f'CHAPTER {n}' in text,(source,n)
 fonts=sorted(set(f for page in reader.pages for f in embedded_fonts(page.get('/Resources'))))
 assert not any(state=='NOT EMBEDDED' for _,state in fonts),(source,fonts)
 records.append({'pdf':str(source.relative_to(ROOT)), 'language':lang,'pages':len(reader.pages),
                 'bytes':source.stat().st_size,'embedded_fonts':fonts})
assert len(records)==12, f'Expected two books and ten chapters; found {len(records)}'
for old,current in json.loads((ROOT/'tools/legacy-downloads.json').read_text()).items():
 assert (ROOT/old).read_bytes()==(ROOT/current).read_bytes(),old
# Both versions and standalone editions should preserve figure numbering.
for lang,book in [('en','IntroImaging'),('fr','IntroImaging-FR')]:
 root=ROOT/'build/latex'/book/'chapters-tex'/lang
 def labels(aux):
  return dict(re.findall(r'\\newlabel\{((?:fig[-:])[^}]+)\}\{\{([^}]+)',aux.read_text()))
 for name in ('01-information','02-images','03-sampling','04-transport','05-neural-networks'):
  assert labels(root/(name+'.aux'))==labels(ROOT/'build/latex'/f'{lang}-{name}'/(f'{lang}-{name}.aux')),(lang,name)
for p in (ROOT/'chapters-tex').rglob('*.tex'):
 s=p.read_text(encoding='utf8');assert s==unicodedata.normalize('NFC',s),p
 assert not re.search('[\u0080-\u009f\ufffd]',s),p
(ROOT/'reports/pdfs.json').write_text(json.dumps(records,indent=2,ensure_ascii=False)+'\n')
print(f'{len(records)} PDFs verified: text, metadata, fonts, dates and numbering.')
