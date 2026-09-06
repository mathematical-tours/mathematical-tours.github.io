from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();mu=np.linspace(1.05,8,18);grid=np.linspace(.4,4,25);scores=np.array([[snr(f,wav_denoise(y,.08*r,'semisoft',mu=m)) for r in grid] for m in mu]);fig,axs=canvas(1,2,height=3.8)
a=axs[0,0].pcolormesh(grid,mu,scores,shading='nearest',cmap='viridis');axs[0,0].set(title='Measured reconstruction SNR',xlabel=r'$T/\sigma$',ylabel=r'$\mu$');fig.colorbar(a,ax=axs[0,0],label='SNR (dB)');curve(axs[0,1],mu,scores.max(axis=1),'Threshold optimized for each μ',r'$\mu$','Best grid SNR (dB)')
finish(fig,__file__,'The heatmap is a genuine two-parameter sweep of semi-soft orthogonal wavelet denoising on one fixed noisy image. The right curve takes a separate maximum over the threshold grid at each mu; it is not a slice through a single globally selected threshold. Tuning uses the clean reference and is explicitly oracle-based.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'mu_grid':mu.tolist(),'threshold_grid':grid.tolist(),'wavelet':'db4','levels':4},arrays={'mu':mu,'threshold_sigma':grid,'snr_db':scores})
