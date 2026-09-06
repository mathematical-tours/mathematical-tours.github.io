from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=.05+.90*flower(128);K=4;rng=np.random.default_rng(20260906);y=f*rng.gamma(K,1/K,f.shape);c=special.digamma(K)-np.log(K);sigma_log=np.sqrt(special.polygamma(1,K));z=np.log(y)-c;sigma_direct=np.sqrt(np.mean(y*y)/(K+1));grid=np.linspace(.5,4,20);direct=[];stab=[]
for r in grid:
    g=stationary_denoise(y,r*sigma_direct,'hard');direct.append((snr(f,g),g,r))
    g=np.exp(stationary_denoise(z,r*sigma_log,'hard'));stab.append((snr(f,g),g,r))
a=max(stab,key=lambda v:v[0]);b=max(direct,key=lambda v:v[0]);fig,axs=canvas(1,3,height=3.6)
for ax,g,title in zip(axs.ravel(),[y,a[1],b[1]],[f'Noisy: {snr(f,y):.1f} dB',f'Log-stabilized: {a[0]:.1f} dB',f'Direct: {b[0]:.1f} dB']):image(ax,g,title)
finish(fig,__file__,'The stabilized target is log(f0): the known mean log multiplier is subtracted before thresholding, and no extra centering constant is added during exponentiation. Both methods use identical noisy data and independently oracle-tuned thresholds. Exponentiating the estimate may introduce bias; the displayed SNR values are actual results, not a claim of unbiasedness or universal superiority.',data_sources=['data/flower.png'],parameters={'seed':20260906,'looks':K,'threshold_grid':grid.tolist(),'clean_intensity':'.05 + .90 flower_luminance'},checks={'stabilized_snr_db':a[0],'direct_snr_db':b[0],'stabilized_threshold_factor':float(a[2]),'direct_threshold_factor':float(b[2])})
