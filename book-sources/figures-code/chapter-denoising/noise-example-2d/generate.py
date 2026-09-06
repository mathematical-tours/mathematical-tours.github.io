from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();fig,axs=canvas(1,2,height=4)
image(axs[0,0],f,r'$f_0$');image(axs[0,1],y,r'$f=f_0+w,\quad w\sim\mathcal{N}(0,\sigma^2 I)$')
finish(fig,__file__,'A single fixed white-Gaussian realization is added to the supplied flower luminance. Both images use the range [0,1]; clipping is for display only. This same clean/noisy pair is reused in the subsequent thresholding comparisons.',data_sources=['data/flower.png'],parameters={'size':256,'sigma':.08,'seed':20260906},checks={'snr_db':snr(f,y)},arrays={'clean':f,'noisy':y})
