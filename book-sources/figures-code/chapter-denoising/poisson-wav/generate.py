from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=flower(128);f/=f.max();peak=30;rng=np.random.default_rng(20260906);counts=rng.poisson(peak*f);y=counts/peak;grid=np.linspace(.5,4,20);sigma=np.sqrt(counts.mean());direct=[];stabilized=[];z=2*np.sqrt(counts+3/8)
for r in grid:
    g=stationary_denoise(counts.astype(float),r*sigma,'hard')/peak;direct.append((snr(f,g),g,r))
    zz=stationary_denoise(z,r,'hard');g=(np.maximum(zz,2*np.sqrt(3/8))**2/4-3/8)/peak;stabilized.append((snr(f,g),g,r))
a=max(direct,key=lambda x:x[0]);b=max(stabilized,key=lambda x:x[0]);fig,axs=canvas(1,3,height=3.6)
for ax,g,title in zip(axs.ravel(),[y,a[1],b[1]],[f'Noisy: {snr(f,y):.1f} dB',f'Direct: {a[0]:.1f} dB',f'Anscombe: {b[0]:.1f} dB']):image(ax,g,title)
finish(fig,__file__,'Both denoisers use the same Poisson counts, stationary wavelet transform, and oracle threshold grid. Direct denoising uses one global count-noise scale sqrt(mean(counts)); stabilized denoising uses Anscombe’s approximate unit scale. The latter projects onto the inverse domain before applying z^2/4-3/8. This direct inverse is not claimed to be unbiased.',data_sources=['data/flower.png'],parameters={'seed':20260906,'maximum_mean_counts':peak,'threshold_grid':grid.tolist()},checks={'direct_snr_db':a[0],'stabilized_snr_db':b[0],'direct_threshold_factor':float(a[2]),'stabilized_threshold':float(b[2])})
