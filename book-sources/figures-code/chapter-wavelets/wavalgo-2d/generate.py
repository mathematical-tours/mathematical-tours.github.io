from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *

f=flower(256)
lo,hi=pywt.dwt(f,'haar',axis=1,mode='periodization');rows=np.concatenate([lo,hi],axis=1)
lo2,hi2=pywt.dwt(rows,'haar',axis=0,mode='periodization');both=np.concatenate([lo2,hi2],axis=0)
fig,axs=canvas(1,3,height=3.7)
image(axs[0,0],f,'Original approximation')
tensor_wavelet_image(axs[0,1],rows,1,'Row analysis',axes=(1,),gain=10,limit=3)
tensor_wavelet_image(axs[0,2],both,1,'Column analysis',gain=10,limit=3)
for x,y,label in [(64,64,'LL'),(192,64,'LH'),(64,192,'HL'),(192,192,'HH')]:
    axs[0,2].text(x,y,label,color='white',ha='center',va='center',fontsize=9,bbox=dict(facecolor=INK,alpha=.7,pad=2,edgecolor='none'))
assert abs(np.sum(f*f)-np.sum(both*both))<1e-8
finish(fig,__file__,'Restores the input photograph, one vertical split after row analysis, and the four quadrants after column analysis. LL/LH/HL/HH describe low/high filtering of the first and second array coordinates. One orthonormal Haar step preserves shape and energy.',parameters={'wavelet':'haar','levels':1},checks={'energy_preserved':True})
