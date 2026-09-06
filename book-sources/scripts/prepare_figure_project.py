#!/usr/bin/env python3
"""Prepare reversible chapter/figure directories; keep published assets unchanged."""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import re
import shutil
from figure_inventory import ROOT, inventory, active, graphics, slug
from tikz_reviews import entries

def display_driver(body):
    return r'''\documentclass[10pt]{book}
\input{book-preamble}
\usepackage[active,tightpage]{preview}
\setlength\PreviewBorder{3pt}
\externaldocument{build/figure-processing/book-external}
\providecommand{\Symb}{s}
\begin{document}
\begin{preview}
\begin{minipage}{\textwidth}\centering
''' + body + r'''
\end{minipage}
\end{preview}
\end{document}
'''

def main():
    if (ROOT/'figure-processing/figure-project.json').exists():
        raise RuntimeError('The figure project is already initialized; use scripts/regenerate_figures.py.')
    figures=inventory(); legacy={e['id']:e for e in entries(ROOT)}
    # Include the four mathematical illustrations without captions and the cover.
    extras=[]
    for source in dict.fromkeys(i['source'] for i in figures):
        text=active((ROOT/source).read_text()); remaining=text
        owned=[i for i in figures if i['source']==source]
        for i in sorted(owned,key=lambda x:x['start'],reverse=True):
            remaining=remaining[:i['start']]+''.join('\n' if c=='\n' else ' ' for c in remaining[i['start']:i['end']])+remaining[i['end']:]
        for a in graphics(remaining):
            ch=owned[0]; name=Path(a['reference']).name
            pos=a['start'];body=text[pos:a['end']]
            extras.append(dict(chapter=ch['chapter'],chapter_title=ch['chapter_title'],chapter_slug=ch['chapter_slug'],
                source=source,start=pos,end=a['end'],line=text.count('\n',0,pos)+1,kind='inline',body=body,
                caption_source='Unnumbered illustration: '+name.replace('-',' '),book_label=None,continued=False,
                assets=graphics(body),context=text[max(0,pos-2500):pos+2000],figure_number=None,book_page=None,
                book_anchor=None,caption_latex='Unnumbered illustration: '+name.replace('-',' '),
                figure_slug=name,id=ch['chapter']+'--'+name))
    cover=ROOT/'figures/wave.png'
    if not cover.exists():cover=next((ROOT/'figures').glob('wave.*'))
    extras.append(dict(chapter='cover',chapter_title='Book Cover',chapter_slug='chapter-book-cover',source='FundationsDataScience.tex',
        start=0,end=0,line=36,kind='cover',body=r'\includegraphics[width=\linewidth]{wave}',
        caption_source='Book cover: wave',book_label=None,continued=False,assets=graphics(r'\includegraphics[width=\linewidth]{wave}'),
        context='A wave motif introduces the signal-processing, optimization, and learning topics of the book.',
        figure_number=None,book_page='Cover',book_anchor=None,caption_latex='Book cover: wave',figure_slug='wave',id='cover--wave'))
    figures.extend(extras)
    folders={};partcounts=defaultdict(int);sources=[]
    for item in figures:
        group=item['figure_number'] or item['id']
        folder=folders.setdefault(group,Path(item['chapter_slug'])/item['figure_slug'])
        partcounts[group]+=1;part=partcounts[group]
        stem='' if part==1 else '-continued'
        code=ROOT/'figures-code'/folder;assetdir=ROOT/'figures'/folder
        code.mkdir(parents=True,exist_ok=True);(assetdir/'original').mkdir(parents=True,exist_ok=True)
        original=item['body']; proposed=item['body'];tikz=[];notes=[];mapped=[]
        for a in item['assets']:
            src=ROOT/a['path'];name=src.name
            dst=assetdir/'original'/name
            if dst.exists() and dst.read_bytes()!=src.read_bytes():
                name=src.stem+'-'+hashlib.sha256(str(src).encode()).hexdigest()[:8]+src.suffix;dst=dst.with_name(name)
            if not dst.exists():shutil.copy2(src,dst)
            replacement=str(dst.relative_to(ROOT/'figures'))
            mapped.append(dict(old=a['reference'],published=replacement,sha256=hashlib.sha256(src.read_bytes()).hexdigest()))
            original=original.replace('{'+a['reference']+'}','{'+replacement+'}')
            if a['command']=='image':
                old=r'\image{'+a['args'][0]+'}{'+a['args'][1]+'}{'+a['args'][2]+'}'
                new=r'\includegraphics[width='+a['args'][1]+r'\linewidth]{'+replacement+'}'
                original=original.replace(old,new)
            if a['reference'].startswith('tikz/'):
                oldtex=src.with_suffix('.tex');newtex=code/oldtex.name
                newtex.write_text(oldtex.read_text().replace('figures/tikz/tikz-preamble.tex','figures-code/tikz-preamble.tex'))
                target=assetdir/(oldtex.stem+'.pdf')
                tikz.append(dict(source=str(newtex.relative_to(ROOT)),output=str(target.relative_to(ROOT))))
                proposed=proposed.replace('{'+a['reference']+'}','{'+str(target.relative_to(ROOT/'figures'))+'}')
                if src.stem in legacy:notes.append(legacy[src.stem]['notes'])
        original_source=code/('original'+stem+'.tex');original_source.write_text(display_driver(original))
        engine='tikz' if tikz and len(tikz)==len(item['assets']) else 'pending'
        # Formula-only Haar figure gets a new plotted construction, not a copy.
        if engine=='tikz':
            (code/('proposed'+stem+'.tex')).write_text(display_driver(proposed))
        item.update(folder=str(folder),part=part,asset_directory=str(assetdir.relative_to(ROOT)),
            code_directory=str(code.relative_to(ROOT)),original_source=str(original_source.relative_to(ROOT)),
            original_pdf=str((assetdir/('original'+stem+'.pdf')).relative_to(ROOT)),
            proposed_pdf=str((assetdir/('proposed'+stem+'.pdf')).relative_to(ROOT)),
            proposed_source=str((code/('proposed'+stem+'.tex')).relative_to(ROOT)) if engine=='tikz' else None,
            engine=engine,tikz_sources=tikz,published_assets=mapped,status='awaiting-generation',
            notes=('Regenerated from the existing editable TikZ construction in its new chapter/figure directory. '+ ' '.join(notes)) if engine=='tikz' else '',
            in_book=False)
        (code/('context'+stem+'.tex')).write_text('% Context from '+item['source']+'; for mathematical reference, not compilation.\n'+item['context'])
    (ROOT/'figures-code/tikz-preamble.tex').write_text((ROOT/'figures/tikz/tikz-preamble.tex').read_text())
    (ROOT/'figure-processing/figure-project.json').write_text(json.dumps(figures,indent=2,ensure_ascii=False)+'\n')
    print(f'Prepared {len(figures)} displays in {len(folders)} figure folders; {sum(i["engine"]=="tikz" for i in figures)} have existing editable constructions.')

if __name__=='__main__':main()
