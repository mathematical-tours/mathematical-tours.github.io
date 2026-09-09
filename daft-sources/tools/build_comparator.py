#!/usr/bin/env python3
"""Build the original/revised figure comparison from the book's actual labels."""
import argparse
import json
from pathlib import Path
import re
import subprocess
from pypdf import PdfReader, PdfWriter, Transformation
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'comparator-figures'


def labels_from_aux(path, seen=None):
    seen=set() if seen is None else seen
    if path in seen:
        return {}
    seen.add(path)
    s=path.read_text()
    labels={k:(num,page) for k,num,page in re.findall(
        r'\\newlabel\{([^}]+)\}\{\{([^}]+)\}\{([^}]+)\}',s)}
    for child in re.findall(r'\\@input\{([^}]+)\}',s):
        candidate=ROOT/child
        if candidate.exists():
            labels.update(labels_from_aux(candidate,seen))
    return labels


def caption_from_source(item):
    text=(ROOT/item['source']).read_text()
    bodies=re.findall(r'\\begin\{figure\}.*?\\end\{figure\}',text,re.S)
    body=next(b for b in bodies if '\\label{'+item['id']+'}' in b)
    label_position=body.index('\\label{'+item['id']+'}')
    start=body.rfind('\\caption{',0,label_position)+len('\\caption{')
    depth=1
    for i in range(start,len(body)):
        if body[i]=='{' and body[i-1]!='\\': depth+=1
        if body[i]=='}' and body[i-1]!='\\': depth-=1
        if depth==0: return body[start:i]
    raise ValueError('Unclosed caption: '+item['id'])


def original_panel(item, legacy):
    """Reassemble original components at their TeX sizes, retaining vector art."""
    cache=ROOT/'build/comparator-inputs'/item['id']
    cache.mkdir(parents=True,exist_ok=True)
    pages=[]
    body=re.sub(r'(?<!\\)%[^\n]*','',item['body'])
    includes=list(re.finditer(r'\\includegraphics\s*(?:\[([^]]*)\])?\s*\{([^}]+)\}',body))
    widths=[]
    forced_breaks=[]
    for i,rel in enumerate(item['graphics']):
        src=legacy/rel
        if src.suffix.lower()=='.eps':
            converted=src.with_name(src.stem+'-eps-converted-to.pdf')
            if converted.exists():
                src=converted
            else:
                pdf=cache/f'{i}.pdf'
                subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-dEPSCrop',
                                '-sDEVICE=pdfwrite',f'-sOutputFile={pdf}',str(src)],check=True)
                src=pdf
        elif not src.suffix:
            src=next((src.with_suffix(ext) for ext in ['.pdf','.png','.jpg','.eps']
                      if src.with_suffix(ext).exists()))
            if src.suffix=='.eps':
                pdf=cache/f'{i}.pdf'
                subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-dEPSCrop',
                                '-sDEVICE=pdfwrite',f'-sOutputFile={pdf}',str(src)],check=True)
                src=pdf
        if src.suffix.lower()!='.pdf':
            pdf=cache/f'{i}.pdf'
            with Image.open(src) as im:
                im.convert('RGB').save(pdf,resolution=150)
            src=pdf
        page=PdfReader(src).pages[0]
        pages.append(page)
        options=includes[i].group(1) or ''
        scale=re.search(r'scale\s*=\s*([.\d]+)',options)
        width=float(page.mediabox.width)*(float(scale[1]) if scale else 1)
        specified=re.search(r'width\s*=\s*([.\d]+)\s*(\\(?:textwidth|linewidth)|cm|mm|in|pt)',options)
        if specified:
            unit=specified[2]
            width=float(specified[1])*({'cm':72/2.54,'mm':72/25.4,'in':72,'pt':72/72.27}.get(unit,468))
        widths.append(width)
        between=body[includes[i-1].end():includes[i].start()] if i else ''
        forced_breaks.append(bool(re.search(r'\\\\|\n\s*\n',between)))
    # The original fullpage layout has a 6.5-inch (468 bp) text block.
    # TeX wraps adjacent graphics when their natural scaled widths do not fit.
    rows=[]; row=[]; used=0; gap=4
    for i,w in enumerate(widths):
        if row and (forced_breaks[i] or used+gap+w>468):
            rows.append(row); row=[]; used=0
        row.append(i); used+=w+(gap if len(row)>1 else 0)
    if row: rows.append(row)
    heights=[w*float(p.mediabox.height)/float(p.mediabox.width) for p,w in zip(pages,widths)]
    row_widths=[sum(widths[i] for i in row)+gap*(len(row)-1) for row in rows]
    row_heights=[max(heights[i] for i in row) for row in rows]
    width=max(row_widths)
    height=sum(row_heights)+8*(len(rows)-1)
    writer=PdfWriter(); dest=writer.add_blank_page(width,height)
    y=height
    for row,rw,rh in zip(rows,row_widths,row_heights):
        y-=rh; x=(width-rw)/2
        for i in row:
            page=pages[i]; scale=widths[i]/float(page.mediabox.width)
            transform=Transformation().translate(-float(page.mediabox.left),-float(page.mediabox.bottom)).scale(scale).translate(x,y)
            dest.merge_transformed_page(page,transform)
            x+=widths[i]+gap
        y-=8
    path=DEST/'before'/(item['id']+'.pdf')
    path.parent.mkdir(parents=True,exist_ok=True)
    writer.write(path)
    return path.relative_to(ROOT).as_posix()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--legacy',type=Path,default=ROOT/'legacy/pre-python-2026-09-08')
    parser.add_argument('--aux',type=Path,default=ROOT/'DAFT-EN.aux')
    parser.add_argument('--partial',action='store_true')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'docs/figure-manifest.json').read_text())
    figures=[f for f in manifest if f['status']=='regenerated']
    if len(figures)!=len(manifest) and not args.partial:
        parser.error('Unfinished figures: use --partial for an explicitly marked working comparison.')
    labels=labels_from_aux(args.aux.resolve())
    DEST.mkdir(exist_ok=True)
    lines=[r'\documentclass[10pt,oneside]{article}',
           r'\input{chapters-tex/shared/book-metadata}',
           r'\usepackage[T1]{fontenc}',r'\usepackage[utf8]{inputenc}',
           r'\usepackage{newpxtext}',r'\usepackage[scaled=.94]{sourcesanspro}',
           r'\input{chapters-tex/shared/daft-macros}',
           r'\let\BookTranspose\transp\let\transp\relax\usepackage{newpxmath}\let\transp\BookTranspose',
           r'\input{ursfs.fd}\DeclareFontShape{U}{rsfs}{m}{n}{<5-7>rsfs5 <7-8>rsfs7 <8->rsfs10}{}',
           r'\usepackage[a4paper,landscape,margin=15mm]{geometry}',
           r'\usepackage{xcolor,hyperref,fancyhdr}',
           r'\definecolor{BookAccent}{HTML}{146C78}',
           r'\hypersetup{colorlinks=true,urlcolor=BookAccent}',
           r'\pagestyle{fancy}\fancyhf{}\fancyfoot[L]{\sffamily\small DAFT --- figure comparison}\fancyfoot[R]{\sffamily\small\thepage}',
           r'\setlength{\parindent}{0pt}\setlength{\emergencystretch}{2em}',
           r'\begin{document}',
           r'{\sffamily\color{BookAccent}\LARGE\bfseries Before and after}\par\medskip',
           r'{\large\BookTitle\par}\smallskip',
           r'\BookAuthor\quad---\quad\BookEditionDate\par\bigskip',
           r'The left panels reproduce the previous English edition. The right panels are regenerated from Python or TikZ and used in the revised book. Figure numbers and captions below are read from the compiled book and its source. The original number is also recorded when needed.\par\bigskip']
    if args.partial:
        lines.append(r'\textbf{Work in progress: '+str(len(figures))+' of '+str(len(manifest))+r' figures regenerated.} The complete comparison will replace this working document.\par')
    audit=[]
    for item in figures:
        num,page=labels[item['id']]
        before=original_panel(item,args.legacy.resolve())
        caption=caption_from_source(item)
        # Resolve references inside captions using the same compiled book labels.
        caption=re.sub(r'\\(eqref|ref)\{([^}]+)\}',
                       lambda m: '('+labels[m[2]][0]+')' if m[1]=='eqref' else labels[m[2]][0],caption)
        lines.extend([r'\clearpage',r'{\sffamily\color{BookAccent}\Large\bfseries Figure '+num+r'}\par\medskip',
                      caption+r'\par\smallskip',
                      r'{\small Book page '+page+r'; original figure '+str(item['original_number'])+r'.}\par\bigskip',
                      r'\begin{minipage}[t]{.485\linewidth}\centering{\sffamily\bfseries Previous edition}\par\medskip',
                      r'\includegraphics[width=\linewidth,height=.66\textheight,keepaspectratio]{'+before+r'}',
                      r'\end{minipage}\hfill\begin{minipage}[t]{.485\linewidth}\centering{\sffamily\bfseries Revised edition}\par\medskip',
                      r'\includegraphics[width=\linewidth,height=.66\textheight,keepaspectratio]{'+item['output']+r'}',
                      r'\end{minipage}'])
        audit.append({'label':item['id'],'number':num,'book_page':page,'original_number':item['original_number'],
                      'caption':caption,'code':item['code'],'output':item['output']})
    lines.append(r'\end{document}')
    (DEST/'comparison.tex').write_text('\n'.join(lines)+'\n')
    (DEST/'index.json').write_text(json.dumps(audit,indent=2)+'\n')
    for _ in range(2):
        result=subprocess.run(['pdflatex','-halt-on-error','-interaction=nonstopmode',
                               '-output-directory=comparator-figures','comparator-figures/comparison.tex'],
                              cwd=ROOT,stdout=subprocess.PIPE,text=True)
        if result.returncode:
            raise RuntimeError(result.stdout[-5000:])
    log=(DEST/'comparison.log').read_text()
    problems=[line for line in log.splitlines() if any(word in line for word in ['Overfull','Underfull','Warning:','Missing character'])]
    print('\n'.join(problems) if problems else 'Comparison compiled without warnings.')
    print(f'{len(figures)} comparisons: {DEST / "comparison.pdf"}')

if __name__=='__main__':
    main()
