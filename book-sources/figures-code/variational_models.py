"""Shared periodic PDE experiments; stopping-time and parameter sweeps."""
from common import *

def tv_flow(y,times,eps=.01):
    targets=np.asarray(times);x=y.copy();time=0.;out=[];dt=.2*eps
    for target in targets:
        while time<target-1e-12:
            step=min(dt,target-time);d=gradient(x);x+=step*divergence(d/np.sqrt(np.sum(d*d,axis=0)+eps**2));time+=step
        out.append(x.copy())
    return out

def flow_experiment():
    f,y=denoising_case(256);heat_times=np.geomspace(.01,5,24);tv_times=np.geomspace(.001,.3,24);sob_params=np.geomspace(.01,10,24);tv_params=np.geomspace(.001,.3,24)
    runs={'Sobolev flow':(heat_times,[heat(y,t) for t in heat_times]),'TV flow':(tv_times,tv_flow(y,tv_times)),
        'Sobolev regularization':(sob_params,[sobolev_denoise(y,p) for p in sob_params]),'TV regularization':(tv_params,[tv_denoise(y,p,400) for p in tv_params])}
    return f,y,runs

def radial_mask(n,K):
    k=fft.fftfreq(n)*n;y,x=np.meshgrid(k,k,indexing='ij');mask=np.zeros((n,n),bool)
    for theta in np.arange(K)*np.pi/K:mask|=abs(-np.sin(theta)*x+np.cos(theta)*y)<=.5
    # Hermitian closure is essential to keep real-valued reconstructions.
    idx=(-np.arange(n))%n;mask|=mask[np.ix_(idx,idx)];return mask

def fourier_tv(y,mask,lam=.02,steps=600):
    x=fft.ifft2(y,norm='ortho').real;bar=x.copy();p=np.zeros((2,*x.shape));tau=sigma=.34
    for _ in range(steps):
        p+=sigma*gradient(bar);p/=np.maximum(1,np.sqrt(np.sum(p*p,axis=0))/lam)
        v=x+tau*divergence(p);new=fft.ifft2((fft.fft2(v,norm='ortho')+tau*y)/(1+tau*mask),norm='ortho').real
        bar=2*new-x;x=new
    return x
