from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


t,f=signal_1d(1024);fig,a=canvas(1,4,height=3.2,width=12)
for ax,M in zip(a.ravel(),[16,32,64,256]):
    ax.plot(t,f,color=MUTED,lw=.8,label='f');ax.plot(t,approx(f,M,'Fourier',True),label=f'M = {M}');ax.set(xlabel='t',title=f'Low-frequency Fourier: M = {M}');ax.grid(True);ax.legend()
finish(fig,__file__,'Low-frequency real Fourier projections of the same piecewise-smooth periodic signal show Gibbs oscillations near jumps. Counts include the constant and individual sine/cosine atoms; ties in frequency are resolved deterministically.',parameters={'N':1024,'M':[16,32,64,256]})
