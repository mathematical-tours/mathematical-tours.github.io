from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from tikz_figure import build
lines=[r'\begin{tikzpicture}[>=Stealth,x=1cm,y=1cm,font=\small,text=Ink,block/.style={draw=Teal,rounded corners=2pt,fill=Teal!5,minimum width=4.1cm,minimum height=.65cm}]']
for k in range(4):
    x=1.2*k
    lines += [rf'\node (e{k}) at ({x},3.1) {{$f[{2*k}]$}};',rf'\node (o{k}) at ({x+5.5},3.1) {{$f[{2*k+1}]$}};']
lines += [r'\node[block] (E) at (1.8,2.25) {$\operatorname{FFT}_4$ (even samples)};',r'\node[block] (O) at (7.3,2.25) {$\operatorname{FFT}_4$ (odd samples)};']
for k in range(4):
    x=1.2*k
    lines += [rf'\draw[->,Teal] (e{k})--({x},2.58);',rf'\draw[->,Teal] (o{k})--({x+5.5},2.58);']
    lines += [rf'\draw[->,Teal] ({x},1.92)--({x},1.45);',rf'\node (E{k}) at ({x},1.2) {{$E_{k}$}};',rf'\draw[->,Teal] ({x+5.5},1.92)--({x+5.5},1.45);',rf'\node (O{k}) at ({x+5.5},1.2) {{$O_{k}$}};']
# Four separate butterflies make each pair of output indices unambiguous.
for k in range(4):
    x=2.45*k+.65
    lines += [rf'\node[draw=Teal,fill=Teal!4,rounded corners,minimum width=2.1cm,minimum height=.72cm] (B{k}) at ({x},-.05) {{$E_{k}\pm\omega_8^{{-{k}}}O_{k}$}};']
    lines += [rf'\draw[->,Teal!70] (E{k}.south)--(B{k}.north west);',rf'\draw[->,Red!65] (O{k}.south)--(B{k}.north east);']
    lines += [rf'\node (p{k}) at ({x-.52},-1.1) {{$\widehat f[{k}]$}};',rf'\node (m{k}) at ({x+.52},-1.1) {{$\widehat f[{k+4}]$}};']
    lines += [rf'\draw[->,Teal] (B{k}.south west)--node[left] {{$+$}} (p{k});',rf'\draw[->,Red] (B{k}.south east)--node[right] {{$-$}} (m{k});']
lines += [r'\end{tikzpicture}']
source=Path(__file__).with_name('figure.tex');source.write_text('\n'.join(lines)+'\n')
build(source,'chapter-discrete-fourier-transform','eight-point-fft')
