from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


t,f,y=noisy_signal();f2,y2=denoising_case();widths=np.linspace(.1,6,60);one=np.array([snr(f,ndimage.gaussian_filter1d(y,s,mode='wrap')) for s in widths]);two=np.array([snr(f2,ndimage.gaussian_filter(y2,s,mode='wrap')) for s in widths]);fig,axs=canvas(1,2,height=3.5)
for ax,values,title in [(axs[0,0],one,'Signal'),(axs[0,1],two,'Image')]:
    curve(ax,widths,values,title,'Gaussian width s (samples/pixels)','SNR (dB)');k=np.argmax(values);ax.plot(widths[k],values[k],'o',color=RED);ax.annotate(f's = {widths[k]:.2f}',(widths[k],values[k]),xytext=(10,-20),textcoords='offset points',fontsize=10)
finish(fig,__file__,'These curves are measured against the known clean signal and flower image for the same observations used in the neighboring examples. The marked grid maximizers are oracle choices for illustration; selecting them requires the clean reference and is not presented as a deployable estimator.',data_sources=['data/flower.png'],parameters={'width_grid':[.1,6,60],'seed':20260906},checks={'signal_optimal_width':float(widths[np.argmax(one)]),'image_optimal_width':float(widths[np.argmax(two)])},arrays={'width':widths,'signal_snr':one,'image_snr':two})
