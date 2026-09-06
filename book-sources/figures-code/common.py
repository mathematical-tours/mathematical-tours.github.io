"""Shared plotting conventions and explicitly defined numerical operators."""
from pathlib import Path
import hashlib
import json
import os
import textwrap
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'build/figure-regeneration/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy import ndimage, fft, signal, special
from scipy.io import wavfile
from PIL import Image
import pywt
from matplotlib.patches import Rectangle

ROOT=Path(__file__).resolve().parents[1]
INK='#203746';TEAL='#146C78';RED='#B94B44';GOLD='#C48A32';BLUE='#4778A3';MUTED='#667780'
COLORS=[TEAL,RED,GOLD,BLUE,'#78578E','#607B43']
DIVERGING=LinearSegmentedColormap.from_list('book-diverging',[TEAL,'#ffffff',RED])
SCALE_COLORS=['#d84b40','#3287ba','#36a06b','#9a61b5','#d99624']
USED_DATA=set()
IMAGE_PREPROCESSING=[]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':12.5,
 'axes.labelsize':11,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,
 'mathtext.fontset':'stix','text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':MUTED,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.spines.top':True,'axes.spines.right':True,
 'axes.prop_cycle':plt.cycler(color=COLORS),'lines.linewidth':1.7,
 'figure.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42,'ps.fonttype':42,
 'axes.titleweight':'normal','axes.titlepad':10,'grid.alpha':.18,'grid.linewidth':.6})

def canvas(rows=1,cols=1,*,height=None,width=9):
    return plt.subplots(rows,cols,figsize=(width,height or 2.8*rows),squeeze=False,layout='constrained')

def image(ax,array,title='',*,signed=False,vmin=None,vmax=None,cmap=None):
    if signed:
        lim=max(float(np.max(np.abs(array))),1e-10);vmin=-lim if vmin is None else vmin;vmax=lim if vmax is None else vmax
    elif array.ndim==2:
        vmin=0 if vmin is None else vmin;vmax=1 if vmax is None else vmax
    obj=ax.imshow(array,cmap=cmap or (DIVERGING if signed else 'gray'),vmin=vmin,vmax=vmax,
                  interpolation='nearest',origin='upper')
    ax.set_title(title);ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_visible(True);spine.set_linewidth(.65)
    return obj

def curve(ax,x,y,title='',xlabel='',ylabel='',label=None,**kwargs):
    ax.plot(x,y,label=label,**kwargs);ax.set(title=title,xlabel=xlabel,ylabel=ylabel);ax.grid(True)
    return ax

def finish(fig,source,notes,*,parameters=None,data_sources=(),checks=None,arrays=None):
    source=Path(source).resolve();folder=source.parent.relative_to(ROOT/'figures-code')
    out=ROOT/'figures'/folder;out.mkdir(parents=True,exist_ok=True)
    for ax in fig.axes:
        title=ax.get_title();width=max(19,int(ax.get_position().width*fig.get_figwidth()*10))
        if '$' not in title:
            ax.set_title('\n'.join(textwrap.fill(line,width=width,break_long_words=False,break_on_hyphens=False) for line in title.split('\n')))
    fig.savefig(out/'proposed.pdf',bbox_inches='tight',pad_inches=.12,metadata={'Title':folder.name.replace('-',' '),'Creator':'Reproducible book figure: '+str(source.relative_to(ROOT))})
    fig.savefig(out/'proposed.png',dpi=150,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
    if getattr(fig,'_book_wavelet_display',False):
        notes+=' Coefficient displays use a shared signed logarithmic grayscale with zero at mid-gray; the coarse approximation is displayed at its image intensity scale. Colored boundaries follow the actual dyadic subband slices, from red at the finest level to blue, green and purple at coarser levels. Display compression does not alter the numerical coefficients.'
    data_sources=sorted(set(map(str,data_sources))|USED_DATA)
    provenance={'source':str(source.relative_to(ROOT)),'notes':notes,'parameters':parameters or {},
        'data_sources':[{'path':str(p),'sha256':hashlib.sha256((ROOT/p).read_bytes()).hexdigest()} for p in data_sources],
        'checks':checks or {},'image_preprocessing':IMAGE_PREPROCESSING}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2,ensure_ascii=False)+'\n')
    if arrays:np.savez_compressed(out/'numerical-data.npz',**arrays)
    else:(out/'numerical-data.npz').unlink(missing_ok=True)
    (source.parent/'README.md').write_text('# '+folder.name.replace('-',' ').title()+'\n\n'+notes+'\n\nRun from the repository root:\n\n```sh\npython '+str(source.relative_to(ROOT))+'\n```\n\nParameters and numerical checks are recorded beside the generated PDF in `provenance.json`.\n')

def flower(n=256,color=False):
    USED_DATA.add('data/flower.png')
    im=Image.open(ROOT/'data/flower.png').convert('RGB')
    entry={'source':'data/flower.png','native_size':list(im.size),'output_size':[n,n],'channels':'RGB' if color else 'Rec. 709 luminance','resampling':'Lanczos'}
    if entry not in IMAGE_PREPROCESSING:IMAGE_PREPROCESSING.append(entry)
    f=np.asarray(im.resize((n,n),Image.Resampling.LANCZOS),dtype=float)/255
    return f if color else np.einsum('ijk,k->ij',f,[.2126,.7152,.0722])

def flower_detail(n=256):
    """Central detail of the supplied flower; crop before resizing."""
    path='data/flower.png';USED_DATA.add(path);im=Image.open(ROOT/path).convert('RGB');w,h=im.size
    bounds=(w//4,h//4,3*w//4,3*h//4);im=im.crop(bounds).resize((n,n),Image.Resampling.LANCZOS)
    entry={'source':path,'crop_pixels':list(bounds),'output_size':[n,n],'grayscale':'Rec. 709 luminance','resampling':'Lanczos'}
    if entry not in IMAGE_PREPROCESSING:IMAGE_PREPROCESSING.append(entry)
    return np.einsum('ijk,k->ij',np.asarray(im,dtype=float)/255,[.2126,.7152,.0722])

PHANTOM_ELLIPSES=np.array([
    [1,.69,.92,0,0,0],[-.8,.6624,.874,0,-.0184,0],
    [-.2,.11,.31,.22,0,-18],[-.2,.16,.41,-.22,0,18],
    [.1,.21,.25,0,.35,0],[.1,.046,.046,0,.1,0],
    [.1,.046,.046,0,-.1,0],[.1,.046,.023,-.08,-.605,0],
    [.1,.023,.023,0,-.606,0],[.1,.023,.046,.06,-.605,0]])

def phantom(n=128):
    """Modified Shepp-Logan's ten analytic ellipses on [-1,1]^2."""
    y,x=np.mgrid[:n,:n]*2/(n-1)-1;y=-y;f=np.zeros((n,n))
    for density,a,b,cx,cy,angle in PHANTOM_ELLIPSES:
        angle=np.deg2rad(angle);u=(x-cx)*np.cos(angle)+(y-cy)*np.sin(angle);v=-(x-cx)*np.sin(angle)+(y-cy)*np.cos(angle)
        f+=density*((u/a)**2+(v/b)**2<=1)
    return f

def phantom_projection(t,theta):
    """Exact continuous Radon projection of the same elliptical phantom."""
    result=np.zeros_like(t,dtype=float)
    for density,a,b,cx,cy,angle in PHANTOM_ELLIPSES:
        alpha=theta-np.deg2rad(angle);radius=np.hypot(a*np.cos(alpha),b*np.sin(alpha));center=cx*np.cos(theta)+cy*np.sin(theta)
        result+=2*density*a*b/radius*np.sqrt(np.maximum(1-((t-center)/radius)**2,0))
    return result

def wavelet_boundaries(ax,slices,*,color_by_scale=True,labels=False):
    """Draw exact packed subband boundaries at pixel edges, finest scale first."""
    for level,detail in enumerate(reversed(slices[1:]),1):
        color=SCALE_COLORS[(level-1)%len(SCALE_COLORS)] if color_by_scale else SCALE_COLORS[0]
        for bounds in detail.values():
            rows,cols=bounds;y0=rows.start or 0;x0=cols.start or 0
            ax.add_patch(Rectangle((x0-.5,y0-.5),cols.stop-x0,rows.stop-y0,fill=False,edgecolor=color,lw=.85,zorder=5))
        if labels:
            r,c=detail['dd'];ax.text(c.stop-3,r.stop-3,f'$D_{{{level}}}$',ha='right',va='bottom',color=color,fontsize=8,
                                   bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=.5),zorder=6)

def wavelet_image(ax,coefficients,slices,title='',*,gain=20,limit=4,labels=False):
    """Gray signed coefficients and an unaltered coarse image, like the book."""
    display=.5+.5*np.sign(coefficients)*np.minimum(np.log1p(gain*np.abs(coefficients))/limit,1)
    coarse=slices[0];a=coefficients[coarse];scale=2**len(slices[1:])
    display[coarse]=np.clip(a/scale,0,1)
    obj=image(ax,display,title,vmin=0,vmax=1)
    wavelet_boundaries(ax,slices,labels=labels)
    ax.figure._book_wavelet_display=True
    return obj

def tensor_boundaries(ax,shape,levels,*,axes=(0,1)):
    """Independent one-dimensional dyadic partitions, not an isotropic pyramid."""
    for level in range(1,levels+1):
        color=SCALE_COLORS[(level-1)%len(SCALE_COLORS)]
        if 1 in axes:ax.axvline(shape[1]/2**level-.5,color=color,lw=.85)
        if 0 in axes:ax.axhline(shape[0]/2**level-.5,color=color,lw=.85)

def coefficient_image(ax,a,title='',*,gain=20,limit=4):
    return image(ax,.5+.5*np.sign(a)*np.minimum(np.log1p(gain*abs(a))/limit,1),title,vmin=0,vmax=1)

def tensor_wavelet_image(ax,a,levels,title='',*,axes=(0,1),gain=20,limit=4):
    display=.5+.5*np.sign(a)*np.minimum(np.log1p(gain*abs(a))/limit,1)
    coarse=tuple(slice(0,a.shape[k]//2**levels if k in axes else a.shape[k]) for k in (0,1))
    display[coarse]=np.clip(a[coarse]/2**(levels*len(axes)/2),0,1)
    image(ax,display,title,vmin=0,vmax=1);tensor_boundaries(ax,a.shape,levels,axes=axes)
    ax.figure._book_wavelet_display=True

def bird():
    fs,y=wavfile.read(ROOT/'data/bird.wav');y=y.astype(float)
    if y.ndim>1:y=y.mean(axis=1)
    return fs,y/max(np.max(np.abs(y)),1)

def signal_1d(n=1024):
    t=np.arange(n)/n
    f=.35+.22*np.sin(2*np.pi*t)+.14*np.sin(10*np.pi*t)+.42*(t>.28)-.55*(t>.62)+.2*(t>.83)
    return t,f

def cartoon(n=256):
    y,x=np.mgrid[:n,:n]/n
    return .16+.15*x+.12*y+.50*((x-.48)**2/.28**2+(y-.51)**2/.34**2<1)+.14*((x>.13)&(x<.3)&(y>.12)&(y<.4))

def snr(clean,estimate):
    return float(20*np.log10(np.linalg.norm(clean)/max(np.linalg.norm(clean-estimate),1e-15)))

def gradient(f):return np.stack((np.roll(f,-1,axis=0)-f,np.roll(f,-1,axis=1)-f))
def divergence(p):return p[0]-np.roll(p[0],1,axis=0)+p[1]-np.roll(p[1],1,axis=1)
def laplacian(f):return divergence(gradient(f))
def tv(f,eps=0):return float(np.sqrt(np.sum(gradient(f)**2,axis=0)+eps**2).sum())

def wavelet_array(f,wavelet='db4',level=4):
    coeff=pywt.wavedecn(f,wavelet,level=level,mode='periodization')
    arr,slices=pywt.coeffs_to_array(coeff);return arr,slices

def wavelet_inverse(a,slices,wavelet='db4'):
    return pywt.waverecn(pywt.array_to_coeffs(a,slices,output_format='wavedecn'),wavelet,mode='periodization')

def hard(a,T):return a*(np.abs(a)>T)
def soft(a,T):return np.sign(a)*np.maximum(np.abs(a)-T,0)
def stein(a,T):return a*np.maximum(1-T*T/np.maximum(a*a,1e-30),0)
def semisoft(a,T,mu):
    if mu==1:return hard(a,T)
    return np.sign(a)*np.minimum(np.abs(a),np.maximum(mu/(mu-1)*(np.abs(a)-T),0))

def wav_denoise(y,T,mode='soft',wavelet='db4',level=4,mu=2,block=1):
    coeff=pywt.wavedecn(y,wavelet,level=level,mode='periodization')
    for detail in coeff[1:]:
        for key,a in detail.items():
            if block>1:
                shape=a.shape;pad=[(0,(-n)%block) for n in shape];ap=np.pad(a,pad,mode='wrap')
                if a.ndim==2:
                    b=ap.reshape(ap.shape[0]//block,block,ap.shape[1]//block,block)
                    energy=np.sum(b*b,axis=(1,3),keepdims=True)
                    rms=np.sqrt(energy/(block*block))
                    factor={'hard':lambda:(rms>T).astype(float),'soft':lambda:np.maximum(1-T/np.maximum(rms,1e-30),0),
                            'stein':lambda:np.maximum(1-T*T/np.maximum(rms*rms,1e-30),0)}[mode]()
                    ap=(b*factor).reshape(ap.shape)
                    detail[key]=ap[:shape[0],:shape[1]]
                else:raise ValueError('Block thresholding is implemented for images.')
            else:detail[key]={'soft':soft,'hard':hard,'stein':stein,'semisoft':lambda a,t:semisoft(a,t,mu)}[mode](a,T)
    return pywt.waverecn(coeff,wavelet,mode='periodization')

def cycle_denoise(y,T,mode='soft',shifts=4,**kwargs):
    result=np.zeros_like(y)
    for i in range(shifts):
        for j in range(shifts):
            v=wav_denoise(np.roll(y,(i,j),(0,1)),T,mode,**kwargs)
            result+=np.roll(v,(-i,-j),(0,1))
    return result/(shifts*shifts)

def stationary_denoise(y,T,mode='hard',wavelet='db4',level=3,mu=2,block=1):
    coeff=pywt.swtn(y,wavelet,level=level,trim_approx=True,norm=False)
    for detail in coeff[1:]:
        for key,a in detail.items():
            if block==1:
                detail[key]={'hard':hard,'soft':soft,'stein':stein,'semisoft':lambda a,t:semisoft(a,t,mu)}[mode](a,T)
            else:
                total=np.zeros_like(a)
                for i in range(block):
                    for j in range(block):
                        v=np.roll(a,(i,j),(0,1));b=v.reshape(v.shape[0]//block,block,v.shape[1]//block,block)
                        rms=np.sqrt(np.mean(b*b,axis=(1,3),keepdims=True))
                        factor={'hard':lambda:(rms>T).astype(float),'soft':lambda:np.maximum(1-T/np.maximum(rms,1e-30),0),
                                'stein':lambda:np.maximum(1-T*T/np.maximum(rms*rms,1e-30),0)}[mode]()
                        total+=np.roll((b*factor).reshape(a.shape),(-i,-j),(0,1))
                detail[key]=total/(block*block)
    return pywt.iswtn(coeff,wavelet,norm=False)

def denoising_case(n=256,sigma=.08,seed=20260906):
    f=flower(n);rng=np.random.default_rng(seed);return f,f+sigma*rng.standard_normal(f.shape)

def noisy_signal(n=1024,sigma=.12,seed=20260906):
    t,f=signal_1d(n);return t,f,f+sigma*np.random.default_rng(seed).standard_normal(n)

def tune_threshold(f,y,sigma,mode='soft',stationary=False,grid=None,**kwargs):
    grid=np.linspace(.25,4.5,24) if grid is None else np.asarray(grid)
    scores=[];best=None
    for r in grid:
        g=(stationary_denoise if stationary else wav_denoise)(y,r*sigma,mode,**kwargs)
        value=snr(f,g);scores.append(value)
        if best is None or value>best[0]:best=(value,g,float(r))
    return best,grid,np.array(scores)

def keep_largest(a,M):
    b=np.zeros_like(a);idx=np.argsort(np.abs(a).ravel(),kind='stable')[-M:]
    b.ravel()[idx]=a.ravel()[idx];return b

def heat(f,t):
    fy=2*np.pi*np.fft.fftfreq(f.shape[0]);fx=2*np.pi*np.fft.fftfreq(f.shape[1])
    lam=4*np.sin(fy[:,None]/2)**2+4*np.sin(fx[None,:]/2)**2
    return fft.ifft2(fft.fft2(f)*np.exp(-t*lam)).real

def sobolev_denoise(y,lam):
    fy=2*np.pi*np.fft.fftfreq(y.shape[0]);fx=2*np.pi*np.fft.fftfreq(y.shape[1])
    eigen=4*np.sin(fy[:,None]/2)**2+4*np.sin(fx[None,:]/2)**2
    return fft.ifft2(fft.fft2(y)/(1+lam*eigen)).real

def tv_denoise(y,lam,steps=250):
    # Chambolle--Pock for .5||x-y||^2 + lam*sum ||Dx||_2; ||D||^2 <= 8.
    x=y.copy();bar=x.copy();p=np.zeros((2,*x.shape));tau=sigma=.34
    for _ in range(steps):
        p+=sigma*gradient(bar);p/=np.maximum(1,np.sqrt(np.sum(p*p,axis=0))/lam)
        new=(x+tau*divergence(p)+tau*y)/(1+tau);bar=2*new-x;x=new
    return x

def inpaint(y,mask,kind='sobolev',steps=500):
    x=y.copy()
    if kind=='sobolev':
        for _ in range(steps):x+=.24*laplacian(x);x[mask]=y[mask]
    elif kind=='tv':
        bar=x.copy();p=np.zeros((2,*x.shape))
        for _ in range(steps):
            p+=.34*gradient(bar);p/=np.maximum(1,np.sqrt(np.sum(p*p,axis=0)))
            new=x+.34*divergence(p);new[mask]=y[mask];bar=2*new-x;x=new
    return x

def huffman_codes(counts):
    import heapq
    heap=[];serial=0
    for symbol,weight in sorted(counts.items()):
        if weight>0:heap.append((float(weight),serial,{symbol:''}));serial+=1
    heapq.heapify(heap)
    while len(heap)>1:
        a,_,ca=heapq.heappop(heap);b,_,cb=heapq.heappop(heap)
        merged={s:'0'+v for s,v in ca.items()}|{s:'1'+v for s,v in cb.items()}
        heapq.heappush(heap,(a+b,serial,merged));serial+=1
    return heap[0][2]

def code_tree(ax,codes,title='',weights=None):
    leaves=sorted(codes.values());prefixes={c[:j] for c in leaves for j in range(len(c)+1)}
    x={p:np.mean([i for i,c in enumerate(leaves) if c.startswith(p)]) for p in prefixes}
    for p in sorted(prefixes,key=len):
        if p:
            ax.plot([x[p[:-1]],x[p]],[-len(p)+1,-len(p)],color=MUTED,lw=1)
            ax.text((x[p[:-1]]+x[p])/2,-len(p)+.55,p[-1],color=MUTED,fontsize=9,
                    bbox={'facecolor':'white','edgecolor':'none','pad':.2},ha='center')
        ax.scatter(x[p],-len(p),s=20 if p not in leaves else 100,color=INK if p not in leaves else TEAL,zorder=3)
    for symbol,code in codes.items():
        label=str(symbol)+'\n'+(code or 'empty')
        if weights:label+='\n'+f'{weights[symbol]:g}'
        ax.text(x[code],-len(code)-.16,label,ha='center',va='top',fontsize=9,color=INK)
    ax.text(x[''],.22,r'$\varnothing$',ha='center');ax.set_title(title)
    ax.set(xlim=(-.5,max(len(leaves)-.5,.5)),ylim=(-max(map(len,leaves))-1.05,.7));ax.set_axis_off()

def diagram_axis(ax):
    ax.set(xlim=(-.035,1.075),ylim=(-.03,1.04));ax.set_axis_off()

def box(ax,x,y,label,width=.14,height=.12,*,color=TEAL,size=12):
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch((x-width/2,y-height/2),width,height,boxstyle='round,pad=.008,rounding_size=.015',
                              edgecolor=color,facecolor='#F1F5F6',linewidth=1.2,zorder=3))
    ax.text(x,y,label,ha='center',va='center',fontsize=size,color=color,zorder=4)

def arrow(ax,start,end,label=None,*,color=MUTED,style='->'):
    ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':style,'color':color,'lw':1.2},zorder=1)
    if label:
        x,y=(np.asarray(start)+end)/2;ax.text(x,y+.025,label,ha='center',va='bottom',fontsize=10,color=color)

def encode_image_at_budget(f,budget,codec):
    import io
    im=Image.fromarray(np.round(np.clip(f,0,1)*255).astype('uint8'))
    choices=[]
    if codec=='JPEG':
        for quality in range(1,96):
            stream=io.BytesIO();im.save(stream,format='JPEG',quality=quality,optimize=True)
            payload=stream.getvalue();choices.append((len(payload),payload,quality))
    elif codec=='JPEG2000':
        low,high=1.,300.
        for _ in range(22):
            ratio=(low+high)/2;stream=io.BytesIO();im.save(stream,format='JPEG2000',irreversible=True,quality_mode='rates',quality_layers=[ratio])
            payload=stream.getvalue();choices.append((len(payload),payload,ratio))
            if len(payload)>budget:low=ratio
            else:high=ratio
    else:raise ValueError(codec)
    eligible=[c for c in choices if c[0]<=budget]
    chosen=max(eligible,key=lambda c:c[0]) if eligible else min(choices,key=lambda c:c[0])
    decoded=np.asarray(Image.open(io.BytesIO(chosen[1])),dtype=float)/255
    return decoded,chosen[1],chosen[2]


def mandrill(n=256):
    path='data/figure-inputs/mandrill.png';USED_DATA.add(path)
    im=Image.open(ROOT/path).convert('RGB').resize((n,n),Image.Resampling.LANCZOS)
    return np.einsum('ijk,k->ij',np.asarray(im,dtype=float)/255,[.2126,.7152,.0722])

def cartoon_frame():
    path='data/figure-inputs/felix-1919.jpg';USED_DATA.add(path)
    return np.asarray(Image.open(ROOT/path).convert('L'),dtype=float)/255
