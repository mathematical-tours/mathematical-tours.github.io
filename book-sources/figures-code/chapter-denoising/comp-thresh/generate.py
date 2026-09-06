from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();fig,axs=canvas(1,3,height=3.7);records={}
methods=[('Stationary hard','hard',True,1),('Orthogonal block Stein','stein',False,4),('Stationary block Stein','stein',True,4)]
for ax,(name,mode,stationary,block) in zip(axs.ravel(),methods):
    best,grid,scores=tune_threshold(f,y,.08,mode,stationary,block=block,level=4);image(ax,best[1],f'{name}\n{best[0]:.1f} dB');records[name]={'threshold_sigma':best[2],'snr_db':best[0]}
T=records['Stationary block Stein']['threshold_sigma']*.08;a=stationary_denoise(y,T,'stein',block=4,level=4);b=np.roll(stationary_denoise(np.roll(y,(1,2),(0,1)),T,'stein',block=4,level=4),(-1,-2),(0,1));err=np.max(abs(a-b));assert err<1e-10
finish(fig,__file__,'All three methods are run and independently oracle-tuned on one observation. Stationary block Stein averages all 4 by 4 block alignments within each undecimated subband, so the combined estimator is equivariant to integer pixel translations. This property is checked numerically; averaging only a few image shifts is not labeled fully translation invariant.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'wavelet':'db4','levels':4,'stationary_block_alignment_average':16},checks={'methods':records,'translation_equivariance_error':float(err)})
