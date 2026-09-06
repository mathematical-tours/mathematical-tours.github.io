from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
t,f=signal_1d(512);cs=pywt.wavedec(f,'haar',level=4,mode='periodization');a=np.concatenate(cs)
fig,axs=canvas(2,3,height=4.8,width=11);axs[0,0].plot(t,f);axs[0,0].set(title='Input signal',xlabel=r'$t$');packed=axs[1,0];packed.plot(a,lw=.7);packed.set(title='Packed Haar coefficients',xlabel='Stored index')
detail_axes=[axs[0,1],axs[1,1],axs[0,2],axs[1,2]]
for ax,c,l in zip(detail_axes,cs[1:],[4,3,2,1]):
    ax.stem(np.arange(len(c)),c,linefmt=TEAL,markerfmt='.',basefmt=' ');ax.set(title=f'Detail level {l}: {len(c)} coefficients',xlabel='Translation index')
    for spine in ax.spines.values():spine.set_color(SCALE_COLORS[l-1]);spine.set_visible(True)
edges=np.r_[0,np.cumsum([len(c) for c in cs])]
for k in range(len(cs)):
    color=SCALE_COLORS[(len(cs)-1-k)%len(SCALE_COLORS)];packed.axvspan(edges[k]-.5,edges[k+1]-.5,color=color,alpha=.10);packed.axvline(edges[k]-.5,color=color,lw=.8)
err=abs(np.sum(a*a)-np.sum(f*f));assert err<1e-10
finish(fig,__file__,'The first column contains the input and complete packed Haar transform. The next two columns enlarge its four detail bands, with matching scale colors. All coefficients come from the same orthonormal transform; coarse coefficients precede details from coarse to fine.',parameters={'wavelet':'haar','samples':512,'levels':4},checks={'parseval_error':float(err)},arrays={'signal':f,'coefficients':a})
