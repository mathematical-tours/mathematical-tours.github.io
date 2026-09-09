from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from related_transform_data import triangle
from hartley import fht
f=triangle();n=f.size;m=1024;nu=np.fft.fftshift(np.fft.fftfreq(m)*n)
F=np.fft.fftshift(np.fft.fft(f,m));H=np.fft.fftshift(fht(np.pad(f,(0,m-n))))
np.testing.assert_allclose(H,F.real-F.imag,atol=1e-12)
fig,axes=plt.subplots(1,2,figsize=(8.7,3.2))
discrete(axes[0],np.arange(n),f);axes[0].plot(f,color=TEAL,lw=.6);axes[0].set(title='Original signal',xlim=(-.3,15.3))
axes[1].plot(nu,H,color=TEAL,label=r'$\mathcal{H}f$')
axes[1].plot(nu,F.real,color=BLUE,ls='--',lw=1.1,label=r'$\operatorname{Re}\widehat f$')
axes[1].plot(nu,F.imag,color=RED,ls=':',lw=.9,label=r'$\operatorname{Im}\widehat f$')
axes[1].set(title='Densely sampled spectra',xlim=(-8,8));axes[1].legend(fontsize=9,loc='upper right')
axes_box(axes[0],'sample index $n$','$f[n]$');axes_box(axes[1],r'frequency $\nu$ (cycles per record)','amplitude')
save(fig,'chapter-related-transforms','hartley-and-fourier')
