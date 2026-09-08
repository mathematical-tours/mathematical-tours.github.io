"""Numerical experiment primitives, shared with the approved Python figures."""
from pathlib import Path
import sys,hashlib,json,functools,itertools,io
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'figures-code'))
import common as C
import numpy as np
from scipy import fft,ndimage,special,signal
from PIL import Image,ImageDraw
import pywt
DATA=ROOT/'myst/live/data'; SEED=20260906
COLORS=['#146c78','#b94b44','#4778a3','#c48a32','#78578e']

def arr(a):
    a=np.asarray(a)
    if not np.all(np.isfinite(a)):raise ValueError('Non-finite numerical output')
    # Significant digits retain tiny objective gaps on logarithmic axes.
    return np.array([float(f'{v:.8g}') for v in a.ravel()]).reshape(a.shape).tolist()

def metric(name,value,unit=''):
    if isinstance(value,(int,float,np.number)):value=f'{value:.3g}'
    return {'label':name,'value':str(value)+(' '+unit if unit else '')}

def state(*panels,metrics=(),note=''):
    return {'panels':list(panels),'metrics':list(metrics),'note':note}

def series(x,y,label='',**kw):return {'x':arr(x),'y':arr(y),'label':label,**kw}

def plot(title,*series_,xlabel='x',ylabel='',xlim=None,ylim=None,logx=False,logy=False,equal=False):
    return {'kind':'plot','title':title,'series':list(series_),'xlabel':xlabel,'ylabel':ylabel,'xlim':xlim,'ylim':ylim,'logx':logx,'logy':logy,'equal':equal}

def table(title,columns,rows):return {'kind':'table','title':title,'columns':columns,'rows':[[str(v) for v in row] for row in rows]}

def image(a,title,*,signed=False,scale=None,boundaries=None,boxes=None):
    a=np.asarray(a)
    if not np.all(np.isfinite(a)):raise ValueError('Non-finite image')
    if signed:
        s=float(scale or np.max(abs(a)) or 1);z=np.clip(a/s,-1,1)
        rgb=np.ones((*a.shape,3))
        pos=np.array([.73,.20,.17]);neg=np.array([.13,.40,.70])
        rgb=np.where((z>=0)[...,None],1-abs(z[...,None])*(1-pos),1-abs(z[...,None])*(1-neg));a=rgb
    pix=np.round(255*np.clip(a,0,1)).astype('uint8')
    identifier=hashlib.sha256(pix.tobytes()+str(pix.shape).encode()).hexdigest()[:20]
    p=DATA/'images'/(identifier+'.webp');p.parent.mkdir(parents=True,exist_ok=True)
    if not p.exists():Image.fromarray(pix).save(p,lossless=True,method=4)
    return {'kind':'image','title':title,'src':'images/'+p.name,'width':pix.shape[1],'height':pix.shape[0], 'boundaries':boundaries or [],'boxes':boxes or []}

def waveimage(a,s,title='Wavelet coefficients'):
    z=.5+.5*np.sign(a)*np.minimum(np.log1p(20*abs(a))/4,1);z[s[0]]=a[s[0]]/2**(len(s)-1)
    boundaries=[]
    for level,d in enumerate(reversed(s[1:])):
        for r,c in d.values():boundaries.append({'x':c.start or 0,'y':r.start or 0,'w':c.stop-(c.start or 0),'h':r.stop-(r.start or 0),'color':C.SCALE_COLORS[level%len(C.SCALE_COLORS)]})
    return image(z,title,boundaries=boundaries)

@functools.lru_cache(maxsize=40)
def input_image(name='Flower',n=128):
    if name=='Flower':return C.flower(n)
    if name=='Mandrill':return C.mandrill(n)
    if name=='Phantom':return C.phantom(n)
    if name=='Cartoon':return C.cartoon(n)
    if name=='Film cartoon':return np.array(Image.fromarray((255*C.cartoon_frame()).astype('uint8')).resize((n,n),Image.Resampling.LANCZOS))/255
    if name=='Smooth':
        y,x=np.mgrid[:n,:n]/n;return .5+.2*np.sin(2*np.pi*x)*np.cos(2*np.pi*y)+.1*np.sin(4*np.pi*y)
    raise KeyError(name)

@functools.lru_cache(maxsize=40)
def input_signal(name='Piecewise smooth',n=256):
    t=np.arange(n)/n
    if name=='Piecewise smooth':return C.signal_1d(n)
    if name=='Smooth':return t,.5+.25*np.sin(2*np.pi*t)+.1*np.cos(8*np.pi*t)
    if name=='Chirp':return t,.5+.35*np.sin(2*np.pi*(3*t+15*t*t))
    if name=='Birdsong':
        fs,y=C.bird();x=np.linspace(0,len(y)-1,n).astype(int);return np.arange(n)/n,y[x]
    raise KeyError(name)

@functools.lru_cache(maxsize=20)
def noise(shape):return np.random.default_rng(SEED).standard_normal(shape)

def slider(id,label,values,default=None):return {'id':id,'label':label,'type':'range','values':list(values),'default':default if default is not None else len(values)//2}
def select(id,label,values,default=0):return {'id':id,'label':label,'type':'select','values':list(values),'default':default}
def images(names=('Flower','Cartoon','Mandrill')):return select('image','Image',names)
def signals():return select('signal','Signal',['Piecewise smooth','Smooth','Chirp'])
def snr(f,g):return min(120,C.snr(f,g))

def wave1(f,w,level):
    c=pywt.wavedec(f,w,mode='periodization',level=level)
    return c,np.concatenate(c)

def wav1_approx(f,w,level):
    c,_=wave1(f,w,level);return pywt.waverec([c[0],*[np.zeros_like(v) for v in c[1:]]],w,mode='periodization')

@functools.lru_cache(maxsize=500)
def transform_image(name,kind,n=128):
    from approximation_models import transform,real_fourier
    f=input_image(name,n)
    if kind=='Fourier':a,inv,r=real_fourier(f);order=np.argsort(r.ravel(),kind='stable')
    elif kind=='DCT':a,inv=transform(f,kind);y,x=np.mgrid[:n,:n];order=np.argsort((x*x+y*y).ravel(),kind='stable')
    else:
        a,s=C.wavelet_array(f,kind,3);inv=lambda b:C.wavelet_inverse(b,s,kind)
        # Coarse-to-fine square ordering, independent of the observed amplitudes.
        y,x=np.mgrid[:n,:n];order=np.argsort(np.maximum(x,y).ravel(),kind='stable')
    assert np.linalg.norm(inv(a)-f)<1e-8
    return a,inv,order

@functools.lru_cache(maxsize=2000)
def approximation_image(name,kind,M,nonlinear=True,n=128):
    a,inv,order=transform_image(name,kind,n);idx=np.argsort(abs(a).ravel(),kind='stable')[-M:] if nonlinear else order[:M]
    b=np.zeros_like(a);b.ravel()[idx]=a.ravel()[idx];g=inv(b)
    assert abs(np.sum((input_image(name,n)-g)**2)-np.sum((a-b)**2))<1e-7
    return g,b
