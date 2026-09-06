from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(1,2,height=3.2)
for ax,wave in zip(axs.ravel(),[False,True]):
    pieces=[(-.25,0,0),(0,.5,1),(.5,1,-1),(1,1.25,0)] if wave else [(-.25,0,0),(0,1,1),(1,1.25,0)]
    for a,b,c in pieces:ax.plot([a,b],[c,c],color=TEAL)
    for a,b,c in pieces[:-1]:ax.plot(b,c,'o',mfc='white',mec=TEAL,ms=5)
    for a,b,c in pieces[1:]:ax.plot(a,c,'o',color=TEAL,ms=5)
    ax.set(title=r'Haar wavelet $\psi$' if wave else r'Haar scaling function $\phi$',xlabel=r'$t$',ylim=(-1.35,1.35),xticks=[0,.5,1],yticks=[-1,0,1]);ax.grid(True)
finish(fig,__file__,'Plots the exact Haar functions at j=0,n=0. Half-open interval conventions are marked with filled and open endpoints, without drawing spurious vertical graph segments at jumps. General atoms follow by translation and multiplication by 2^(-j/2) at scale 2^j; both displayed functions have unit L2 norm.',checks={'scaling_norm_squared':1,'wavelet_norm_squared':1,'wavelet_integral':0})
