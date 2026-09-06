#!/usr/bin/env python3
"""Build the pending figure audit, retaining the author-review numbering history."""
from pathlib import Path
import argparse
import json
import shutil
from figure_tex import tex_text, book_labels, braced
from build_book import compile_document

ROOT=Path(__file__).resolve().parents[1]
BUILD=ROOT/'build/figure-regeneration'

PREAMBLE=r'''\documentclass[10pt,openany]{book}
\input{book-preamble}
\geometry{a3paper,landscape,left=16mm,right=16mm,top=14mm,bottom=14mm}
\externaldocument{build/figure-regeneration/comparator-external}[../FundationsDataScience.pdf]
\hypersetup{pdftitle={\BookTitle: New Figure Reconstructions},filecolor=BookAccent}
\renewcommand{\eqdef}{\ensuremath{\stackrel{\mbox{\normalfont\tiny def.}}{=}}}
\newsavebox{\AuditHeading}\newsavebox{\AuditNotes}\newlength{\AuditHeight}
\newcommand{\AuditFigure}[9]{%
 \begin{lrbox}{\AuditHeading}\begin{minipage}{\linewidth}
 {\sffamily\small\color{BookMuted}\BookTitle\hfill RECONSTRUCTION AUDIT\par}\vspace{3mm}
 {\sffamily\LARGE\bfseries\color{BookAccent}#2\hfill\normalsize\mdseries #3\par}\vspace{2mm}
 {\large #4\par}\vspace{2mm}
 {\small\sffamily\color{BookMuted}#8\par}
 \end{minipage}\end{lrbox}%
 \begin{lrbox}{\AuditNotes}\begin{minipage}{\linewidth}
 \hrule\vspace{2mm}{\small\textbf{Mathematical and generation notes.} #7\par}
 \vspace{2mm}{\scriptsize\sffamily\color{BookMuted}#9\par}
 \end{minipage}\end{lrbox}%
 \setlength{\AuditHeight}{\dimexpr\textheight-\ht\AuditHeading-\dp\AuditHeading-\ht\AuditNotes-\dp\AuditNotes-22mm\relax}%
 \noindent\usebox{\AuditHeading}\par\vspace{3mm}
 \noindent\begin{minipage}[t]{.485\linewidth}
 {\large\sffamily\bfseries\color{BookMuted}Current book\par}\vspace{2mm}
 \parbox[c][\AuditHeight][c]{\linewidth}{\centering\includegraphics[width=\linewidth,height=\AuditHeight,keepaspectratio]{#5}}
 \end{minipage}\hfill\begin{minipage}[t]{.485\linewidth}
 {\large\sffamily\bfseries\color{BookAccent}Proposed reconstruction\par}\vspace{2mm}
 \parbox[c][\AuditHeight][c]{\linewidth}{\centering\includegraphics[width=\linewidth,height=\AuditHeight,keepaspectratio]{#6}}
 \end{minipage}\par\vspace{3mm}\noindent\usebox{\AuditNotes}\par\vfill
 {\scriptsize\sffamily\color{BookMuted}Audit identifier: \texttt{\detokenize{#1}}\hfill\thepage\par}\clearpage}
\begin{document}\pagestyle{empty}\setlength{\parindent}{0pt}
\begin{titlepage}\sffamily\color{BookInk}\vspace*{12mm}
{\Large\BookTitle\par}\vspace{15mm}
{\fontsize{38}{44}\selectfont\bfseries New Figure Reconstructions\par}\vspace{6mm}
{\color{BookAccent}\rule{40mm}{1.2pt}\par}\vspace{10mm}
{\Large\BookAuthor\par}\vspace{3mm}{\BookAffiliation\par}\vspace{5mm}{\BookEditionDate\par}\vfill
\begin{minipage}{.8\linewidth}\large
Each page places a current book figure beside its new reproducible reconstruction.
The figure numbers and complete captions identify the accompanying main book.
Continued figures keep their original numbers and identify the continued display.
Unnumbered illustrations and the cover are listed separately.
Accepted reconstructions have been integrated into the reading editions and are omitted here.
This volume contains only figures for which further revisions were requested.
Where numbering changed, the previous audit number appears beside the current book number.\medskip\par
These revised reconstructions remain proposals for review. The left column shows their current book illustrations.
Use a figure number, together with a panel or continuation when needed, to request changes.
The footer identifies the corresponding generation-code directory.\medskip\par
The notes distinguish measured data, controlled numerical experiments, and schematic constructions.
Parameters, random seeds, source data, and numerical checks are stored beside each figure.
\end{minipage}\vspace{12mm}\end{titlepage}
'''

def build_comparisons(preview=False,chapter_filter=None):
    BUILD.mkdir(parents=True,exist_ok=True)
    from build_book import chapter_list
    labels=book_labels(ROOT,ROOT/'build',chapter_list())
    items=json.loads((ROOT/'figure-processing/figure-project.json').read_text())
    # Refresh exact numbers, pages and captions from the stabilized reading edition.
    for item in items:
        if item.get('review_state') in ('removed','merged'):continue
        if item['book_label']:
            info=labels[item['book_label']];item.update(figure_number=info[0],book_page=info[1],caption_latex=info[2],book_anchor=info[3])
        elif item['figure_number']:
            aux=(ROOT/'build/book/chapters'/f"{item['chapter']}.aux").read_text()
            for line in aux.splitlines():
                if not line.startswith(r'\@writefile{lof}{\contentsline {figure}'):continue
                body,_=braced(line,len(r'\@writefile{lof}'));pos=body.index(r'\numberline')+len(r'\numberline');num,pos=braced(body,pos);cap,pos=braced(body,pos);page,pos=braced(body,pos+1);anchor,_=braced(body,pos)
                if num==item['figure_number']:item.update(book_page=page,book_anchor=anchor,caption_latex=item['caption_source']);break
    # Keep the unnumbered illustrations with their owning chapters. The initial
    # inventory appended these after the numbered figures; that also duplicated
    # chapter bookmark destinations when the audit returned to earlier chapters.
    chapter_order={'cover':0, **{name:n for n,name in enumerate(chapter_list(),1)}}
    items.sort(key=lambda item: chapter_order[item['chapter']])
    (ROOT/'figure-processing/figure-project.json').write_text(json.dumps(items,indent=2,ensure_ascii=False)+'\n')
    items=[i for i in items if i.get('comparison_visible',i['engine']!='tikz')]
    if chapter_filter:items=[i for i in items if i['chapter'] in chapter_filter.split(',')]
    ready=[i for i in items if (ROOT/i['proposed_pdf']).exists()]
    if len(ready)!=len(items) and not preview:raise RuntimeError(f'{len(items)-len(ready)} figures still need reconstruction.')
    # Tell xr-hyper that these are already five-field hyperref labels. Without
    # the marker it appends compatibility fields and corrupts \hyperref URLs.
    external=r'\HyperFirstAtBeginDocument{}'+'\n'+(ROOT/'build/figure-processing/book-external.aux').read_text();pages=[];index=[];chapter=None
    for item in ready:
        label='audit-original-'+item['id']
        if item['book_anchor']:
            external+=r'\newlabel{'+label+'}{{'+item['figure_number']+'}{'+item['book_page']+'}{'+item['caption_latex']+'}{'+item['book_anchor']+'}{}}'+'\n'
        if item['chapter']!=chapter:
            chapter=item['chapter'];pages.append(r'\pdfbookmark[0]{'+item['chapter_title']+'}{chapter-'+chapter+'}')
        heading=('Figure '+item['figure_number']+(' (continued)' if item['continued'] else '')) if item['figure_number'] else 'Unnumbered illustration'
        pages.append(r'\pdfbookmark[1]{'+heading+'}{audit-'+item['id']+'}')
        if item['book_anchor']:
            heading=r'\hyperref['+label+']{'+heading+'}'
        page='Main book, page '+item['book_page'] if item['book_page'] else 'Unnumbered: '+item['chapter_title']
        previous=item.get('review_number')
        if previous and previous!=item['figure_number']:page+=' (previous Figure '+previous+')'
        if item['id']=='denoising--filtering-optimal-1d':page+='; previous Figures 7.9–7.10'
        notes=item['notes'] or 'Generation notes are recorded in the figure source.'
        args=[item['id'],heading,page,item['caption_latex'],item['original_pdf'],item['proposed_pdf'],tex_text(notes),
              item['chapter_title'],r'Code: \texttt{\detokenize{'+item['code_directory']+'}}']
        pages.append('\\AuditFigure'+'\n'.join('{'+x+'}' for x in args))
        index.append({k:item[k] for k in ['id','chapter','figure_number','book_page','book_label','book_anchor','caption_latex','part','original_pdf','proposed_pdf','code_directory','in_book','review_number','review_state']})
        index[-1]['comparison_pdf_page']=len(index)+1
    (BUILD/'comparator-external.aux').write_text(external)
    preamble=PREAMBLE
    if not items:
        preamble=preamble.replace('New Figure Reconstructions','Figure Review Complete')
        preamble=preamble.replace(r'\begin{minipage}{.8\linewidth}\large',r'\begin{minipage}{.8\linewidth}\large\raggedright')
        start=preamble.index('Each page places a current book figure')
        end=preamble.index(r'\end{minipage}',start)
        accepted=sum(i.get('review_state')=='accepted' for i in json.loads((ROOT/'figure-processing/figure-project.json').read_text()))
        preamble=preamble[:start]+(
            f'All {accepted} accepted reconstructions are now included in the main book and its independent chapters. '
            'There are no pending figure comparisons.\\medskip\\par\n'
            'The original cover and the illustrations explicitly retained by the author remain unchanged. '
            'The previous comparison volume and its figure index are preserved in '
            r'\texttt{archive/2026-09-06-reviewed-figures/}.\medskip\par'+'\n'
            'The complete decision ledger and editable figure sources retain the original audit identifiers and current book numbering.\n'
        )+preamble[end:]
    (BUILD/'complete-comparisons.tex').write_text(preamble+'\n'.join(pages)+'\n\\end{document}\n')
    report=compile_document('complete-figure-comparisons',BUILD/'complete-comparisons.tex',BUILD/'comparison-build')
    (BUILD/'comparison-report.json').write_text(json.dumps(report,indent=2)+'\n')
    if report['diagnostics']:raise RuntimeError('\n'.join(report['diagnostics']))
    target=BUILD/'preview-comparisons.pdf' if preview else ROOT/'figure-processing/figure-comparisons.pdf'
    shutil.copy2(report['pdf'],target)
    (BUILD/'complete-figure-index.json').write_text(json.dumps(index,indent=2,ensure_ascii=False)+'\n')
    if not preview:shutil.copy2(BUILD/'complete-figure-index.json',ROOT/'figure-processing/figure-index.json')
    print(f'{len(ready)} comparisons, {report["pages"]} pages: {target}')

    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--preview',action='store_true');p.add_argument('--chapter');a=p.parse_args()
    build_comparisons(a.preview,a.chapter)

if __name__=='__main__':main()
