from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();fig,axs=canvas(1,2,height=3.8);records={}
for mode in ['hard','soft','stein']:
    best,grid,scores=tune_threshold(f,y,.08,mode,block=4);curve(axs[0,0],grid,scores,'Rules on 4 × 4 blocks',r'$T/\sigma$','SNR (dB)',label=mode);records[mode]=best[0]
sizes=[1,2,4,8];maxima=[]
for s in sizes:maxima.append(tune_threshold(f,y,.08,'stein',block=s)[0][0])
curve(axs[0,1],sizes,maxima,'Stein rule: block-size comparison','Block side length','Best grid SNR (dB)',marker='o');axs[0,1].set_xticks(sizes);axs[0,0].legend()
finish(fig,__file__,'The left panel applies hard, soft, and Stein factors to the block RMS, with the block-size normalization included. The right panel independently optimizes T/sigma for each block side length. Blocks remain inside one scale and orientation, and every comparison uses the same noisy image.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'block_sizes':sizes},checks={'best_snr_per_rule':records,'best_snr_per_size':maxima})
