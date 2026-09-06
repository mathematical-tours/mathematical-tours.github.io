from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


t,f=signal_1d(256);fig,axs=canvas(1,3,height=3.1,width=11)
for ax,level in zip(axs.ravel(),[1,2,3]):
    cs=pywt.wavedec(f,'haar',level=level,mode='periodization');arr=np.concatenate(cs);ax.plot(arr,lw=.8)
    pos=0
    for k,c in enumerate(cs):
        ax.axvspan(pos,pos+len(c),color=COLORS[k%len(COLORS)],alpha=.09);ax.text(pos+len(c)/2,ax.get_ylim()[1]*.86,'A' if k==0 else f'D{level-k+1}',ha='center',fontsize=10);pos+=len(c)
    ax.set(title=f'{level} decomposition '+('step' if level==1 else 'steps'),xlabel='Storage index');assert arr.size==f.size;assert abs(np.sum(arr**2)-np.sum(f**2))<1e-10
finish(fig,__file__,'The storage arrays are recomputed after one, two and three Haar steps. The active approximation occupies the leftmost block and is replaced by a smaller approximation followed by its detail. Earlier details stay in place. Length and energy are preserved at every step.',parameters={'samples':256,'wavelet':'haar','levels':[1,2,3]},checks={'length_preserved':True,'energy_preserved':True})
