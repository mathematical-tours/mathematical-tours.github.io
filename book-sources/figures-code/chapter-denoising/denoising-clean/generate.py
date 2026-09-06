from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();g=stationary_denoise(y,2.7*.08,'hard');fig,axs=canvas(1,3,height=3.5)
for ax,v,title in zip(axs.ravel(),[f,y,g],['Clean image',f'Noisy: {snr(f,y):.1f} dB',f'Denoised: {snr(f,g):.1f} dB']):image(ax,v,title)
finish(fig,__file__,'The estimate is computed by hard thresholding stationary Daubechies-4 wavelet details at 2.7 sigma and applying the exact inverse stationary transform. The coarsest approximation is retained. Measured SNR is shown for the noisy and reconstructed images.',data_sources=['data/flower.png'],parameters={'sigma':.08,'seed':20260906,'threshold_sigma':2.7,'wavelet':'db4','stationary_levels':3},checks={'noisy_snr_db':snr(f,y),'denoised_snr_db':snr(f,g)})
