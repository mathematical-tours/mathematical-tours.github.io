from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
f=flower(256);bounds=(72,56,152,136);x0,y0,x1,y1=bounds
crop=f[y0:y1,x0:x1];small=np.asarray(Image.fromarray((255*crop).astype('uint8')).resize((20,20),Image.Resampling.BOX))/255
fs,y=bird();window=max(16,int(.01*fs));energy=np.convolve(y*y,np.ones(window)/window,mode='same');center=int(np.argmax(energy));lo=max(0,center-int(.001*fs));hi=min(len(y),center+int(.001*fs))
fig,axs=canvas(1,4,height=3.3,width=12)
image(axs[0,0],f,'Flower image');axs[0,0].add_patch(Rectangle((x0-.5,y0-.5),x1-x0,y1-y0,fill=False,edgecolor=RED,lw=2))
image(axs[0,1],small,'Magnified region: 20 × 20 pixels')
for sp in axs[0,1].spines.values():sp.set_color(RED);sp.set_linewidth(2)
t=np.arange(len(y))/fs;axs[0,2].plot(t,y,lw=.55);axs[0,2].axvspan(lo/fs,hi/fs,color=RED,alpha=.25);axs[0,2].set(title='Recorded birdsong',xlabel='Time (s)',ylabel='Amplitude')
axs[0,3].plot(1000*(t[lo:hi]-t[center]),y[lo:hi],lw=1);axs[0,3].set(title='Short-time detail',xlabel='Time from excerpt center (ms)',ylabel='Amplitude')
for sp in axs[0,3].spines.values():sp.set_color(RED)
finish(fig,__file__,'The red spatial box selects precisely the flower crop shown as 20 by 20 area-averaged pixels. The sound panels use the supplied bird recording and a 2 ms excerpt around a strong chirp; the red interval marks this exact excerpt. No temporal subsampling is applied.',data_sources=['data/bird.wav'],parameters={'crop_pixels':bounds,'zoom_grid':[20,20],'sample_rate':int(fs),'excerpt_samples':[lo,hi]})
