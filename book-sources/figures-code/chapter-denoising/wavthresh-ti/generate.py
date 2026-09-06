from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();fig,axs=canvas(height=3.6,width=7);records={}
for mode in ['hard','soft','stein']:
    best,grid,values=tune_threshold(f,y,.08,mode,True,level=4);curve(axs[0,0],grid,values,'Stationary wavelet thresholding',r'$T/\sigma$','SNR (dB)',label=mode);records[mode]={'threshold_sigma':best[2],'snr_db':best[0]}
axs[0,0].legend();finish(fig,__file__,'Every point is measured after thresholding the stationary Daubechies-4 coefficients and reconstructing the image. Hard, soft, and Stein rules share the same observation, levels, and threshold grid. The redundant transform retains all translations; no old SNR curves are reused.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'wavelet':'db4','levels':4},checks=records)
