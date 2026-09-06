from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fs,y=bird();t=np.arange(len(y))/fs;f=flower(256)
fig,axs=canvas(1,3,height=3.3)
curve(axs[0,0],t,y,'Sound: one time coordinate',r'$t$ (seconds)','Amplitude')
image(axs[0,1],f,'Image: two spatial coordinates')
ax=axs[0,2];ax.set_axis_off();ax.set_title('Video: space and time')
for j in range(4):
    frame=np.roll(f,8*j,axis=1)
    ax.imshow(frame,cmap='gray',vmin=0,vmax=1,extent=(j*.12,1+j*.12,j*.15,1+j*.15),zorder=4-j)
ax.set(xlim=(-.04,1.6),ylim=(-.22,1.55));ax.annotate('time',xy=(1.4,-.13),xytext=(.15,-.13),arrowprops={'arrowstyle':'->','color':TEAL},color=TEAL)
finish(fig,__file__,'The waveform uses the supplied bird recording with its actual sample rate and duration. The grayscale image is the supplied flower image. The video is an explicitly simulated sequence of translated flower frames; it illustrates d=3, not an additional measured recording.',data_sources=['data/bird.wav','data/flower.png'],parameters={'sample_rate':int(fs),'video_frames':4,'translation_pixels_per_frame':8},checks={'audio_duration_seconds':len(y)/fs})
