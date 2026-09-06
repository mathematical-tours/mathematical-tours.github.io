from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


t,f=signal_1d();fig,a=canvas(1,4,height=3.3,width=12)
for ax,M in zip(a.ravel(),[16,32,64,128]):ax.plot(t,f,color=MUTED,lw=.8);ax.plot(t,approx(f,M),label=f'M={M}');ax.set(title=f'Largest db4 coefficients: M = {M}',xlabel='t');ax.grid(True)
finish(fig,__file__,'Nonlinear orthonormal wavelet approximations retain the largest 16,32,64,128 coefficients of the same 1024-sample signal. Fine details survive primarily near jumps.',parameters={'N':1024,'M':[16,32,64,128],'wavelet':'db4','levels':4})
