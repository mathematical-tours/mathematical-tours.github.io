from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();fig,axs=canvas(1,2,height=4);records={}
for ax,mode,stationary in [(axs[0,0],'soft',False),(axs[0,1],'hard',True)]:
    best,grid,values=tune_threshold(f,y,.08,mode,stationary,level=4);name='Stationary hard' if stationary else 'Orthogonal soft';image(ax,best[1],f'{name}: {best[0]:.1f} dB');records[name]={'threshold_sigma':best[2],'snr_db':best[0]}
finish(fig,__file__,'Compares orthogonal soft thresholding with stationary hard thresholding on the same noisy flower. Both use Daubechies-4 filters, four levels, and the same oracle threshold grid. The stationary transform uses every translation and is inverted exactly; it is not an average over only a few shifts.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'wavelet':'db4','levels':4},checks=records)
