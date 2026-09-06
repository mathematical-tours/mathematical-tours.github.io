from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *

from curvelet_frame import windows,analyze,synthesize,denoise
f=flower_detail(128);sigma=.08;y=f+sigma*np.random.default_rng(20260906).standard_normal(f.shape)
bank,info=windows(128);c=analyze(f,bank);err=np.linalg.norm(synthesize(c,bank)-f)/np.linalg.norm(f)
energy=sum(np.sum(abs(v)**2) for v in c);assert err<1e-12;assert abs(energy/np.sum(f*f)-1)<1e-12
grid=np.linspace(1,4,16);runs=[(snr(f,g:=denoise(y,sigma,r,bank)),g,r) for r in grid]
best=max(runs,key=lambda z:z[0]);wb,_,_=tune_threshold(f,y,sigma,'hard',True,grid=grid);fig,axs=canvas(1,3,height=3.8)
for ax,g,title in zip(axs.ravel(),[y,wb[1],best[1]],[f'Noisy: {snr(f,y):.1f} dB',f'Translation-invariant wavelets: {wb[0]:.1f} dB',f'Curvelets: {best[0]:.1f} dB']):image(ax,g,title)
finish(fig,__file__,'Retains the original close-up/noisy/wavelet/curvelet comparison using a detail of the supplied flower. Both denoisers are recomputed from one seeded noisy observation, with oracle thresholds selected on the same grid and SNR measured against the known clean flower detail. Curvelets use the documented periodic oversampled Parseval frame, with smooth parabolic frequency windows and all pixel translations; this is not the fast wrapping implementation.',parameters={'size':128,'sigma':sigma,'seed':20260906,'curvelet_windows':len(bank),'threshold_grid':grid.tolist()},checks={'frame_inverse_relative_error':float(err),'parseval_ratio':float(energy/np.sum(f*f)),'curvelet_threshold_factor':float(best[2]),'wavelet_threshold_factor':wb[2],'curvelet_snr_db':best[0],'wavelet_snr_db':wb[0]})
