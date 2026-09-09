from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from chirp import zoom_spectrum,scaled_fourier
n=128;beta=2.23;t=np.arange(n);f=np.cos(2*np.pi*beta*t/n)+.2*(np.random.default_rng(52).random(n)-.5)
nu=np.linspace(1.7,2.8,1001);zoom=zoom_spectrum(f,nu[0],nu[1]-nu[0],nu.size)
estimate=nu[np.argmax(abs(zoom))];b=int(np.floor(estimate));L=n*b/estimate;r=round(L)
adjusted=scaled_fourier(f[:r],r/L)
fig,axes=plt.subplots(2,2,figsize=(9,6));a,bax,c,d=axes.ravel()
a.plot(t,f,color=TEAL);a.set(title='(a) Observed signal',xlim=(0,n-1));axes_box(a,'sample index $n$','amplitude')
k=np.arange(32);bax.plot(k,abs(np.fft.fft(f)[:32])/n,'o-',ms=3,color=TEAL);bax.set(title='(b) Coarse spectrum',xlim=(0,31));axes_box(bax,'cycles per original record','magnitude / $N$')
c.plot(nu,abs(zoom)/n,color=TEAL);c.axvline(beta,color=RED,ls=':',label=rf'true $\beta={beta}$')
c.axvline(estimate,color=BLUE,ls='--',label=rf'peak $\widehat\beta={estimate:.3f}$');c.set(title='(c) Spectral zoom',xlim=(nu[0],nu[-1]));c.legend(fontsize=8,loc='lower left');axes_box(c,'cycles per original record','magnitude / $N$')
d.plot(k,abs(adjusted[:32])/r,'o-',ms=3,color=TEAL);d.set(title=rf'(d) Adjusted record: $r={r}$',xlim=(0,31));axes_box(d,'adjusted frequency index $k$','magnitude / $r$')
for ax in [bax,d]:ax.set_ylim(-.012,.58)
save(fig,'chapter-related-transforms','spectral-zoom-and-adjustment')
