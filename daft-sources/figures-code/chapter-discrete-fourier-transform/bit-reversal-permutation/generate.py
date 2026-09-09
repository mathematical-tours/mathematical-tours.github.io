from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from tikz_figure import build
from fft_algorithms import bit_reversal_indices
lines=[r'\begin{tikzpicture}[>=Stealth,x=1cm,y=.52cm,font=\small,text=Ink]']
permutation=bit_reversal_indices(8)
lines += [r'\node[above] at (0,.5) {input index};',r'\node[above] at (5,.5) {output position};']
for i,j in enumerate(permutation):
    lines += [rf'\node[anchor=east] (a{i}) at (0,{-i}) {{$f[{i}]$\quad\texttt{{{i:03b}}}}};',rf'\node[anchor=west] (b{i}) at (5,{-i}) {{\texttt{{{i:03b}}}\quad$f[{j}]$}};']
for i,j in enumerate(permutation):
    color='Ink!30' if i==j else ('Teal' if i<j else 'Red')
    lines += [rf'\draw[->,{color},line width=.8pt] (a{i}.east)--(b{j}.west);']
lines += [r'\end{tikzpicture}']
source=Path(__file__).with_name('figure.tex');source.write_text('\n'.join(lines)+'\n')
build(source,'chapter-discrete-fourier-transform','bit-reversal-permutation')
