from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(1,2,height=4)
code_tree(axs[0,0],{format(i,'03b'):format(i,'03b') for i in range(8)},'Complete depth-3 binary tree')
codes={0:'001',1:'01',2:'1',3:'000'};code_tree(axs[0,1],codes,'Prefix code on four symbols')
assert all(not b.startswith(a) for a in codes.values() for b in codes.values() if a!=b)
finish(fig,__file__,'The full tree contains all eight length-three words. The pruned tree uses exactly the original code 0→001, 1→01, 2→1, 3→000. Edge labels append bits; the leaves show the source symbol and its codeword. The prefix property is checked.',parameters={'code':codes},checks={'prefix_free':True,'kraft_sum':sum(2**(-len(c)) for c in codes.values())})
