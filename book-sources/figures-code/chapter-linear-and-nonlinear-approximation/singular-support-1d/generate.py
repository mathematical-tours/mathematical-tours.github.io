from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


fig,a=canvas(1,2,height=3.7);t,f=signal_1d(512);cs=pywt.wavedec(f,'haar',level=6,mode='periodization');ax=a[0,0];ax.plot(t,f);ax.set(title='Piecewise-smooth signal',xlabel='t');ax=a[0,1]
for level in range(1,7):
 c=cs[-level];width=2**level/512;left=np.arange(len(c))*width;center=left+width/2;sing=np.zeros(len(c),bool)
 for jump in [.28,.62,.83,0]:sing|=(left<=jump)&(jump<left+width)
 ax.scatter(center,np.full(len(c),level),c=np.where(sing,RED,TEAL),s=8+60*abs(c)/max(abs(c).max(),1e-10))
ax.set(title='Haar supports meeting a jump (red)',xlabel='Support center',ylabel='Level (coarser upward)');finish(fig,__file__,'Every dot represents an actual Haar detail coefficient. Color is determined by whether its dyadic support intersects a known jump, and area by magnitude. Periodicity adds the boundary jump. Haar is used so support intersection is exact and visible.',parameters={'N':512,'levels':6,'wavelet':'haar'})
