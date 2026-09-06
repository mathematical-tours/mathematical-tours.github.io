from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


from collections import Counter
n=48;values=np.clip(np.floor(1.5+1.4*np.sin(np.linspace(0,3*np.pi,n))+.5),0,3).astype(int);d=np.diff(values);counts=Counter(d.tolist());codes=huffman_codes(counts)
fig,axs=canvas(1,4,height=3.4,width=13)
axs[0,0].step(np.arange(n),values,where='mid',color=TEAL);axs[0,0].plot(np.arange(n),values,'.',color=TEAL);axs[0,0].set(title='Original values',xlabel='Sample index',ylabel=r'$f_i$',yticks=range(4))
axs[0,1].stem(np.arange(1,n),d,linefmt=TEAL,markerfmt='o',basefmt=' ');axs[0,1].set(title='Successive differences',xlabel='Sample index',ylabel=r'$d_i=f_i-f_{i-1}$',yticks=[-1,0,1])
x=np.arange(-1,4);raw=Counter(values.tolist());axs[0,2].bar(x-.17,[raw[k]/n for k in x],width=.34,label='values',color=TEAL);axs[0,2].bar(x+.17,[counts[k]/len(d) for k in x],width=.34,label='differences',color=RED);axs[0,2].set(title='Empirical symbol frequencies',xlabel='Symbol',ylabel='Probability',xticks=x);axs[0,2].legend()
code_tree(axs[0,3],codes,'Prefix tree for the differences')
reconstructed=np.r_[values[0],values[0]+np.cumsum(d)];assert np.array_equal(reconstructed,values)
finish(fig,__file__,'A fixed discrete signal is differenced, and both histograms and the Huffman tree are recomputed from its actual counts. Reconstruction retains the initial value and cumulatively sums the differences. Difference symbols are more concentrated at zero for this slowly varying example; this is an example, not a universal decrease of marginal entropy.',parameters={'samples':n,'initial_value':int(values[0]),'difference_code':{str(k):v for k,v in codes.items()}},checks={'exact_reconstruction':True},arrays={'values':values,'differences':d})
