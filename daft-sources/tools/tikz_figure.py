"""Compile a figure's own TikZ source to a standalone vector PDF."""
from pathlib import Path
import shutil
import subprocess
import pypdfium2 as pdfium

ROOT = Path(__file__).resolve().parents[1]


def build(source, chapter, name):
    source = Path(source).resolve()
    work = ROOT / 'build/tikz' / chapter / name
    work.mkdir(parents=True, exist_ok=True)
    wrapper = work / 'figure.tex'
    wrapper.write_text(r'''\documentclass[border=5pt]{standalone}
\usepackage[T1]{fontenc}
\usepackage{newpxtext,newpxmath}
\usepackage{tikz}
\usetikzlibrary{arrows.meta,calc,positioning}
\definecolor{Ink}{HTML}{203746}
\definecolor{Teal}{HTML}{146C78}
\definecolor{Red}{HTML}{C34E42}
\definecolor{Blue}{HTML}{3478AB}
\begin{document}
''' + source.read_text() + '\n'+r'\end{document}'+'\n')
    command = shutil.which('pdflatex') or '/Library/TeX/texbin/pdflatex'
    process = subprocess.run([command, '-halt-on-error', '-interaction=nonstopmode',
                              'figure.tex'], cwd=work, text=True, capture_output=True)
    log = (work / 'figure.log').read_text(errors='replace')
    if process.returncode or any(s in log for s in ['Overfull', 'Underfull', 'LaTeX Warning']):
        raise RuntimeError(log[-8000:])
    folder = ROOT / 'figures' / chapter / name
    folder.mkdir(parents=True, exist_ok=True)
    shutil.copy2(work/'figure.pdf', folder/'figure.pdf')
    preview = ROOT / 'build/figure-previews' / chapter
    preview.mkdir(parents=True, exist_ok=True)
    doc = pdfium.PdfDocument(str(folder/'figure.pdf'))
    page = doc[0]
    page.render(scale=2).to_pil().save(preview/(name+'.png'))
    page.close(); doc.close()
