"""Synthesis sparse reconstruction with verified orthogonal/Parseval operators."""
from common import *

def wavelet_operator(shape,stationary=False):
    zero=np.zeros(shape)
    if not stationary:
        _,s=wavelet_array(zero,'db2',3)
        return lambda f:wavelet_array(f,'db2',3)[0],lambda a:wavelet_inverse(a,s,'db2')
    template=pywt.swtn(zero,'db2',level=3,trim_approx=True,norm=True);keys=[sorted(d) for d in template[1:]]
    def analysis(f):
        c=pywt.swtn(f,'db2',level=3,trim_approx=True,norm=True);return np.stack([c[0]]+[d[k] for d,kk in zip(c[1:],keys) for k in kk])
    def synthesis(a):
        c=[a[0]];pos=1
        for kk in keys:c.append({k:a[pos+j] for j,k in enumerate(kk)});pos+=len(kk)
        return pywt.iswtn(c,'db2',norm=True)
    return analysis,synthesis

def ista(y,A,At,lam,L=1,steps=500,accelerated=False,x0=None,history=False):
    x=np.zeros_like(At(y)) if x0 is None else x0.copy();z=x.copy();t=1.;values=[];iterates=[]
    for k in range(steps):
        new=soft(z-At(A(z)-y)/L,lam/L)
        if accelerated:tnew=(1+np.sqrt(1+4*t*t))/2;z=new+(t-1)/tnew*(new-x);t=tnew
        else:z=new
        x=new
        if history:values.append(.5*np.linalg.norm(A(x)-y)**2+lam*np.sum(abs(x)));iterates.append(x.copy())
    return (x,np.array(values),np.array(iterates)) if history else x

def blur_symbol(shape,width=1.2):
    grids=np.meshgrid(*[fft.fftfreq(n) for n in shape],indexing='ij');return np.exp(-2*np.pi**2*width**2*sum(k*k for k in grids))
def convolution(x,H):return fft.ifftn(fft.fftn(x)*H).real

def spikes_case():
    n=256;k=fft.fftfreq(n)*n;h=(1-(k/3)**2)*np.exp(-.5*(k/3)**2);h-=h.mean();H=fft.fft(h);H/=abs(H).max();h=fft.ifft(H).real;x=np.zeros(n);x[[27,65,116,155,201]]=[1,-.8,.7,1,-.6];y=convolution(x,H)+.002*np.random.default_rng(20260906).standard_normal(n);return x,h,H,y

def deblur_case(n=256):
    f=flower_detail(n);H=blur_symbol(f.shape);y=convolution(f,H)+.015*np.random.default_rng(20260906).standard_normal(f.shape);return f,H,y

def sparse_deblur(y,H,lam,stationary=False,steps=400):
    W,S=wavelet_operator(y.shape,stationary);A=lambda c:convolution(S(c),H);At=lambda v:W(convolution(v,np.conj(H)));delta=np.zeros_like(y);delta.flat[0]=1
    weights=np.linalg.norm(W(delta),axis=tuple(range(1,W(delta).ndim)),keepdims=True) if stationary else 1
    return S(ista(y,A,At,lam*weights,steps=steps,accelerated=True))

def deblur_sweep():
    f,H,y=deblur_case();grid=np.geomspace(1e-5,1,31);k=fft.fftfreq(f.shape[0]);e=4*np.sin(np.pi*k[:,None])**2+4*np.sin(np.pi*k[None,:])**2;Y=fft.fft2(y)
    runs={'Squared norm':[fft.ifft2(H*Y/(H*H+l)).real for l in grid],
          'Sobolev':[fft.ifft2(H*Y/(H*H+l*e)).real for l in grid],
          'Orthogonal wavelets':[sparse_deblur(y,H,l) for l in grid]}
    return f,H,y,grid,runs

def lq_plot(ax,q):
    t=np.linspace(0,2*np.pi,801);x=np.cos(t);y=np.sin(t);r=(abs(x)**q+abs(y)**q)**(-1/q) if np.isfinite(q) else 1/np.maximum(abs(x),abs(y));ax.fill(r*x,r*y,color=TEAL,alpha=.16);ax.plot(r*x,r*y,color=TEAL);ax.axhline(0,color=MUTED,lw=.6);ax.axvline(0,color=MUTED,lw=.6);ax.set(aspect='equal',xlim=(-1.25,1.25),ylim=(-1.25,1.25),title=rf'$q={q:g}$' if np.isfinite(q) else r'$q=\infty$');ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1])
