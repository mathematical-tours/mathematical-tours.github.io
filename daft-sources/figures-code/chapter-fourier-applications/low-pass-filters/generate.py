from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from image_processing import smooth_lowpass_kernel,filter_response
n=65;k=np.arange(-(n//2),n//2+1);frequency=np.fft.fftshift(np.fft.fftfreq(n))
fig,axes=plt.subplots(3,3,figsize=(9,6.3))
for row,epsilon in enumerate([0,.4,1]):
    g=smooth_lowpass_kernel(n,epsilon);H=np.fft.fftshift(np.fft.fft(g)).real
    discrete(axes[row,0],k,np.fft.fftshift(g));axes[row,0].set_ylabel(rf'$\epsilon={epsilon:g}$')
    axes[row,1].plot(2*frequency,H,'o',ms=2.5,color=TEAL)
    omega,curve=filter_response(g,2048);assert np.max(abs(curve.imag))<1e-13
    axes[row,2].plot(omega/np.pi,curve.real,color=TEAL)
    for col,ax in enumerate(axes[row]):axes_box(ax,'kernel index' if col==0 else r'$\omega/\pi$')
    axes[row,0].set(xlim=(-20,20),ylim=(-.14,.56));axes[row,1].set(xlim=(-1,1),ylim=(-.15,1.15));axes[row,2].set(xlim=(-1,1),ylim=(-.15,1.15))
for ax,title in zip(axes[0],['Impulse response','Prescribed samples','FIR frequency response']):ax.set_title(title)
save(fig,'chapter-fourier-applications','low-pass-filters')
