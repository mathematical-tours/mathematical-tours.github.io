"""Reproducible approximation models, orthogonal transforms and adaptive meshes."""
from common import *
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

def models(n=256):
    y,x=np.mgrid[:n,:n]/n
    smooth=.4+.18*np.sin(2*np.pi*x)*np.cos(2*np.pi*y)+.13*np.cos(4*np.pi*x+2*np.pi*y)
    return {'Smooth':smooth,'Cartoon':cartoon(n),'Flower':flower(n),'Mandrill':mandrill(n)}

def real_fourier(f):
    # tensor product of the real, orthonormal trigonometric basis; one real coefficient per atom.
    n=f.shape[0];t=np.arange(n)/n;columns=[np.ones(n)/np.sqrt(n)];freq=[0]
    for k in range(1,(n+1)//2):
        columns.extend([np.sqrt(2/n)*np.cos(2*np.pi*k*t),np.sqrt(2/n)*np.sin(2*np.pi*k*t)]);freq.extend([k,k])
    if n%2==0:columns.append((-1.)**np.arange(n)/np.sqrt(n));freq.append(n//2)
    B=np.array(columns).T
    a=B.T@f if f.ndim==1 else B.T@f@B
    inverse=(lambda a:B@a) if f.ndim==1 else (lambda a:B@a@B.T)
    r=np.array(freq) if f.ndim==1 else np.hypot(np.array(freq)[:,None],np.array(freq)[None,:])
    return a,inverse,r

def transform(f,kind='Wavelet'):
    if kind=='Fourier':return real_fourier(f)[:2]
    if kind=='DCT':return fft.dctn(f,norm='ortho'),lambda a:fft.idctn(a,norm='ortho')
    if kind=='Local DCT':
        n=f.shape[0];bs=16
        def run(a,inverse=False):
            blocks=a.reshape(n//bs,bs,n//bs,bs);fn=fft.idctn if inverse else fft.dctn
            return fn(blocks,axes=(1,3),norm='ortho').reshape(n,n)
        return run(f),lambda a:run(a,True)
    a,s=wavelet_array(f,'db4',4);return a,lambda b:wavelet_inverse(b,s,'db4')

def approx(f,M,kind='Wavelet',linear=False):
    if kind=='Fourier':
        a,inv,r=real_fourier(f);b=keep_largest(a,M) if not linear else np.where(np.argsort(np.argsort(r.ravel(),kind='stable'),kind='stable').reshape(a.shape)<M,a,0)
    else:
        a,inv=transform(f,kind)
        if linear:
            # packed coarse coefficients occupy the upper-left square (2-D) or first interval (1-D).
            b=np.zeros_like(a)
            if a.ndim==1:b[:M]=a[:M]
            else:
                side=int(np.sqrt(M));assert side*side==M;b[:side,:side]=a[:side,:side]
        else:b=keep_largest(a,M)
    return inv(b)

def edge_mesh(n=256,ring=70):
    t=np.linspace(0,2*np.pi,ring,endpoint=False);pts=[]
    for d in [-.006,.006]:pts.extend(np.c_[.48+(.28+d)*np.cos(t),.51+(.34+d)*np.sin(t)])
    x,y=np.meshgrid(np.linspace(0,1,9),np.linspace(0,1,9));pts.extend(np.c_[x.ravel(),y.ravel()]);points=np.array(pts);return points,Delaunay(points)

def greedy_mesh(f,count=240):
    n=f.shape[0];yy,xx=np.mgrid[:n,:n];grid=np.c_[xx.ravel(),yy.ravel()];pts=[(0,0),(0,n-1),(n-1,0),(n-1,n-1)];chosen=np.zeros_like(f,dtype=bool)
    for x,y in pts:chosen[y,x]=True
    for _ in range(count-4):
        points=np.array(pts);vals=f[points[:,1],points[:,0]];tri=Delaunay(points);g=LinearNDInterpolator(tri,vals)(grid).reshape(f.shape)
        err=(f-g)**2;err[chosen]=-1;y,x=np.unravel_index(np.nanargmax(err),f.shape);pts.append((x,y));chosen[y,x]=True
    points=np.array(pts);tri=Delaunay(points);g=LinearNDInterpolator(tri,f[points[:,1],points[:,0]])(grid).reshape(f.shape)
    return points,tri,g
