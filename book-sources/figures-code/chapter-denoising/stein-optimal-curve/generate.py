from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();fig,axs=canvas(height=3.7,width=7);records={}
for mode in ['stein','hard','soft']:
    best,grid,scores=tune_threshold(f,y,.08,mode,True,level=4);curve(axs[0,0],grid,scores,'Stein and comparison shrinkage rules',r'$T/\sigma$','SNR (dB)',label=mode);records[mode]={'threshold_sigma':best[2],'snr_db':best[0]}
axs[0,0].legend();finish(fig,__file__,'Stationary Daubechies-4 thresholding is evaluated with Stein, hard, and soft nonlinearities on identical data. SNR values and grid maxima are recomputed, so the figure records the result of this experiment rather than asserting that one method is universally best.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'levels':4},checks=records)
