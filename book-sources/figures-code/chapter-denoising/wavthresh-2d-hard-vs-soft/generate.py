from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();fig,axs=canvas(height=3.7,width=7);records={}
for mode in ['hard','soft']:
    best,grid,values=tune_threshold(f,y,.08,mode);curve(axs[0,0],grid,values,'Orthogonal wavelet thresholding',r'$T/\sigma$','SNR (dB)',label=mode);records[mode]={'threshold_sigma':best[2],'snr_db':best[0]}
axs[0,0].legend();finish(fig,__file__,'Hard and soft thresholding use the same noisy flower, orthonormal Daubechies-4 basis, and threshold grid. Each SNR value is measured from an actual reconstructed image. The maxima are reference-dependent experimental optima, not general rankings of the shrinkage rules.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'wavelet':'db4','levels':4},checks=records)
