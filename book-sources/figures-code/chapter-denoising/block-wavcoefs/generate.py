from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();best,grid,scores=tune_threshold(f,y,.08,'stein',block=4);a,s=wavelet_array(y,'db4',4);b,_=wavelet_array(best[1],'db4',4);fig,axs=canvas(1,3,height=3.6)
wavelet_image(axs[0,0],a,s,'Noisy wavelet coefficients');wavelet_image(axs[0,1],b,s,'After 4 × 4 block shrinkage');image(axs[0,2],best[1],f'Reconstruction: {best[0]:.1f} dB')
finish(fig,__file__,'Disjoint 4 by 4 blocks are formed separately within every orthogonal detail subband. Each block shares the Stein factor max(1-T^2/E_B,0), where E_B is its mean squared coefficient magnitude. The threshold is oracle-selected on the common grid, and the displayed coefficient arrays come from the same reconstruction.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'block_side':4,'threshold_sigma':best[2],'wavelet':'db4','levels':4},checks={'snr_db':best[0]})
