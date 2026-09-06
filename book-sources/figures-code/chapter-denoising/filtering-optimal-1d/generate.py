from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
t,f,y=noisy_signal();f2,y2=denoising_case();widths=np.linspace(.1,6,60)
s=widths[np.argmax([snr(f,ndimage.gaussian_filter1d(y,v,mode='wrap')) for v in widths])];s2=widths[np.argmax([snr(f2,ndimage.gaussian_filter(y2,v,mode='wrap')) for v in widths])]
g=ndimage.gaussian_filter1d(y,s,mode='wrap');g2=ndimage.gaussian_filter(y2,s2,mode='wrap');fig=plt.figure(figsize=(13,3.1),layout='constrained');gs=fig.add_gridspec(1,4,width_ratios=[1.3,1.3,1,1])
for col,v,title in [(0,y,f'Noisy signal: {snr(f,y):.1f} dB'),(1,g,f'Filtered: {snr(f,g):.1f} dB')]:
    ax=fig.add_subplot(gs[0,col]);ax.plot(t,f,color=MUTED,lw=.6);curve(ax,t,v,title,r'$t$');ax.set_ylim(-.2,1.4)
image(fig.add_subplot(gs[0,2]),y2,f'Noisy image: {snr(f2,y2):.1f} dB');image(fig.add_subplot(gs[0,3]),g2,f'Filtered: {snr(f2,g2):.1f} dB')
finish(fig,__file__,'Combines the previously separate signal and image comparisons on one row. Both oracle filter widths use exactly the observations and width grid of the preceding SNR plot. Clean references are used for selection and SNR measurement; the image panels share one fixed intensity range.',parameters={'seed':20260906,'signal_width':float(s),'image_width':float(s2)},checks={'signal_snr_db':snr(f,g),'image_snr_db':snr(f2,g2)})
