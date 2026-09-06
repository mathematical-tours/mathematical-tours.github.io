from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
t,f=signal_1d(1024);fig,axs=canvas(3,2,height=6.4,width=10);records={}
for col,w in enumerate(['haar','db4']):
    curve(axs[0,col],t,f,'Input for '+('Haar' if w=='haar' else 'Daubechies 4'),r'$t$');cs=pywt.wavedec(f,w,level=4,mode='periodization')
    for row,L in enumerate([1,3],1):
        c=cs[-L];centers=(np.arange(len(c))+.5)*2**L/len(f)
        axs[row,col].stem(centers,c,linefmt=TEAL,markerfmt='.',basefmt=' ')
        for jump in [.28,.62,.83]:axs[row,col].axvline(jump,ls='--',color=RED,lw=.8)
        axs[row,col].set(title=f'Detail level {L}',xlabel='Coefficient location on the periodic grid',xlim=(0,1));records[w+str(L)]=c
for row in [1,2]:
    lim=max(max(abs(ax.get_ylim()[0]),abs(ax.get_ylim()[1])) for ax in axs[row]);
    for ax in axs[row]:ax.set_ylim(-lim,lim)
finish(fig,__file__,'Haar and Daubechies-4 coefficients of the same piecewise-smooth signal are compared at two common detail levels. Dashed lines identify the input jumps; the coefficient grid uses the same periodic packing convention. Daubechies support spans more samples and is not identified with the one-cell Haar support. Rows share vertical ranges.',parameters={'samples':1024,'wavelets':['haar','db4'],'detail_levels':[1,3]},arrays=records)
