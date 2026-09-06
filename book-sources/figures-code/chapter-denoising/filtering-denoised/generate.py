from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
t,f,y=noisy_signal();n=len(f);d=np.minimum(np.arange(n),n-np.arange(n));s=3.;h=np.exp(-d*d/(2*s*s));h/=h.sum();Y=fft.fft(y);H=fft.fft(h);g=fft.ifft(Y*H).real;G=fft.fft(g);freq=fft.fftshift(fft.fftfreq(n))
fig,axs=canvas(2,3,height=4.8,width=11)
for col,v,V,title in [(0,y,Y,'Noisy signal'),(1,fft.fftshift(h),H,'Gaussian kernel'),(2,g,G,'Filtered signal')]:
    curve(axs[0,col],t if col!=1 else np.arange(-n//2,n//2),v,title,'Time' if col!=1 else 'Sample offset')
    if col==1:axs[0,col].set_xlim(-20,20)
    curve(axs[1,col],freq,np.log(np.maximum(abs(fft.fftshift(V)),1e-10)),title+' spectrum','Cycles per sample',r'$\log|\widehat f|$' if col!=1 else r'$\log|\widehat h|$')
err=np.max(abs(G-Y*H));assert err<1e-10
finish(fig,__file__,'The top row contains signal, kernel and filtered signal; the bottom row contains their corresponding Fourier spectra. All six panels belong to one cyclic convolution g=y*h, G=YH. The Gaussian has unit discrete mass; the plotting floor changes only logarithmic display.',parameters={'samples':n,'sigma_kernel_samples':s,'seed':20260906},checks={'convolution_fourier_residual':float(err),'kernel_mass':float(h.sum())})
