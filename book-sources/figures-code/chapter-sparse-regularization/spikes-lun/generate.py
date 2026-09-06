from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


f,h,H,y=spikes_case();A=lambda x:convolution(x,H);At=lambda x:convolution(x,H.conj());grid=np.geomspace(.0002,.015,16);runs=[(snr(f,g:=ista(y,A,At,l,steps=1800,accelerated=True)),g,l) for l in grid];best=max(runs,key=lambda v:v[0]);inv=np.zeros_like(H);ok=abs(H)>1e-12;inv[ok]=1/H[ok];pinv=convolution(y,inv);fig,a=canvas(2,3,height=5,width=11);n=np.arange(len(f));panels=[(fft.fftshift(h),'Transmitted second-Gaussian-derivative pulse'),(fft.fftshift(abs(H)),'Magnitude of the Fourier symbol'),(f,'Sparse reflectivity'),(y,'Noisy convolution'),(pinv,'Discrete pseudoinverse (note vertical scale)'),(best[1],f'Sparse reconstruction: {best[0]:.1f} dB')]
for ax,(z,title) in zip(a.ravel(),panels):curve(ax,n,z,title,'Index','Amplitude')
finish(fig,__file__,'All six panels arise from one discrete periodic spike experiment. The pseudoinverse discards only numerically zero Fourier modes (tolerance 1e-12), so its noise amplification is visible. Sparse recovery uses the true adjoint and an oracle-selected lambda on a stated grid.',parameters={'N':256,'seed':20260906,'noise_sigma':.002,'lambda_grid':grid.tolist(),'selected_lambda':float(best[2]),'pseudoinverse_tolerance':1e-12},checks={'selected_snr_db':best[0]})
