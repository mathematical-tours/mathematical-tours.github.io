from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


t,f=signal_1d(1024);fig,axs=canvas(3,1,height=6.2)
for ax,s in zip(axs.ravel(),[1,6,20]):
    g=ndimage.gaussian_filter1d(f,s,mode='wrap');ax.plot(t,f,color=MUTED,lw=.65,alpha=.65,label='original');curve(ax,t,g,rf'Gaussian width $\sigma={s}/1024$',r'$t$','Amplitude',label='filtered');ax.legend(loc='upper right',ncols=2)
finish(fig,__file__,'The three smoothed curves are computed from one piecewise smooth signal by convolution with normalized discrete Gaussian kernels, using periodic boundary conditions. The displayed standard deviations are in physical time units. The original signal is retained as a thin guide so widening the filter visibly trades detail for smoothness.',parameters={'samples':1024,'sigma_in_samples':[1,6,20],'boundary':'periodic'})
