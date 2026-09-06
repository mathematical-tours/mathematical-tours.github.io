from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();fig,axs=canvas(1,2,height=4);records={}
for ax,mode in zip(axs.ravel(),['hard','soft']):
    best,grid,values=tune_threshold(f,y,.08,mode);image(ax,best[1],f'{mode.title()}: T/σ = {best[2]:.2f}\nSNR {best[0]:.1f} dB');records[mode]={'threshold_sigma':best[2],'snr_db':best[0]}
finish(fig,__file__,'Shows the actual reconstructions at the hard/soft SNR maxima from the preceding plot, with each method tuned on the same grid and clean reference. Both displays use identical intensity bounds.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'wavelet':'db4'},checks=records)
