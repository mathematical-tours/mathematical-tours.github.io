from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from tikz_figure import build
nodes=[]; edges=[]
for weight in range(5):
    row=[i for i in range(16) if i.bit_count()==weight]
    for j,i in enumerate(row):
        x=(j-(len(row)-1)/2)*1.55
        nodes.append(r"\node[vertex] (v"+str(i)+") at ("+str(x)+","+str(weight*1.15)+") {"+format(i,'04b')+"};")
for i in range(16):
    for bit in range(4):
        j=i^(1<<bit)
        if i<j: edges.append(r"\draw[Teal!55,line width=.65pt] (v"+str(i)+") -- (v"+str(j)+");")
# Nodes are redrawn after the edges so that no line crosses a label.
preamble=r"\begin{tikzpicture}[vertex/.style={draw=Teal,fill=white,rounded corners=5pt,inner xsep=5pt,inner ysep=3pt,font=\small\ttfamily,text=Ink}]"
source=Path(__file__).with_name('figure.tex')
source.write_text(preamble+'\n'+'\n'.join(nodes+edges+nodes)+r"\end{tikzpicture}")
build(source,'chapter-duality','boolean-cube')
