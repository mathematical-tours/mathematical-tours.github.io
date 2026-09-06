from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


u=np.linspace(-4,4,1801);fig,axs=canvas(1,2,height=3.8)
for mu in [1.25,2,4]:curve(axs[0,0],u,semisoft(u,1,mu),'Semi-soft thresholding',r'$x/T$',r'$S_T(x)/T$',label=rf'$\mu={mu}$')
axs[0,0].plot(u,u,'--',color=MUTED,lw=.8);axs[0,0].legend();curve(axs[0,1],u,stein(u,1),'Stein thresholding',r'$x/T$',r'$S_T(x)/T$',label='Stein');axs[0,1].plot(u,soft(u,1),'--',label='soft',color=RED);axs[0,1].plot(u,u,':',color=MUTED,label='identity');axs[0,1].legend()
finish(fig,__file__,'Plots the exact semi-soft rules with transition interval [T,mu T] and the exact Stein attenuation max(1-T^2/x^2,0)x, defining the origin to be zero. Soft thresholding is shown for comparison: its large-amplitude shrinkage tends to T, whereas Stein shrinkage tends to zero.',parameters={'mu_values':[1.25,2,4],'T_units':1},checks={'stein_at_zero':float(stein(np.array([0.]),1)[0]),'semisoft_at_threshold':float(semisoft(np.array([1.]),1,2)[0])})
