"""Noise, regularization, inverse problems, and sparse recovery experiments."""
from numerics import *
import common as C
from variational_models import tv_flow,radial_mask,fourier_tv
from sparse_models import wavelet_operator,ista,blur_symbol,convolution,sparse_deblur


def controls(family,view):
    if family=='noise':return [slider('sigma','Noise standard deviation σ',[0,.02,.04,.08,.12,.2]),images(('Flower','Cartoon'))]
    if family=='gaussian':return [slider('width','Gaussian width (pixels)',[0,.25,.5,.75,1,1.5,2,3,5]),images(('Flower','Cartoon'))]
    if family=='threshold':
        if view in ['semisoft','semisoft-optimal']:return [slider('factor','Threshold T / σ',[.5,1,1.5,2,2.5,3,4]),slider('mu','Transition μ',[1.1,1.5,2,3,5])]
        if 'block' in view or view=='comp-thresh':return [slider('factor','Threshold T / σ',[.5,1,1.5,2,2.5,3,4]),select('block','Block side',[1,2,4])]
        if view=='cycle-spin-principle':return [slider('shift','Translation (pixels)',[0,1,2,3,4,5,6,7]),select('method','Method',['Orthogonal','Cycle spinning'])]
        return [slider('factor','Threshold T / σ',[.5,1,1.5,2,2.5,3,4]),select('method','Threshold',['Soft','Hard','Stein']),images(('Flower','Cartoon'))]
    if family=='poisson':return [slider('photons','Peak photon count',[1,2,4,8,16,32,64,128]),images(('Flower','Cartoon'))]
    if family=='multiplicative':return [slider('looks','Independent looks K',[1,2,4,8,16,32]),images(('Flower','Cartoon'))]
    if family=='gradient' and view in ['grad-discrete','discrete-operators']:return [slider('gain','Display contrast gain',[.5,1,2,4,8,16]),images()]
    if family=='gradient':return [slider('eps','Smoothing ε',[.001,.003,.01,.03,.1,.3]),images()]
    if family=='pde':return [slider('strength','Evolution / regularization strength',[0,.01,.03,.06,.1,.2,.4]),select('method','Model',['Heat flow','TV flow','Sobolev regularization','TV regularization']),images(('Flower','Cartoon'))]
    if family=='inpainting':return [slider('steps','Iterations',[0,10,30,80,200,500]),select('method','Prior',['Sobolev','TV','Wavelet']),select('mask','Missing region',['Square','Stripes','Random'])]
    if family=='tomography':
        if view=='tomo-principle':return [slider('angle','Projection angle θ (degrees)',[0,15,30,45,60,75,90,120,150,180])]
        return [slider('lines','Radial Fourier lines',[4,8,13,20,32,48,64]),select('method','Reconstruction',['Pseudoinverse','TV'])]
    if family=='sparse':
        if view=='sparsity-wav':return [slider('factor','Display / selection threshold',[0,.02,.05,.1,.2,.4]),images()]
        return [slider('lam','Regularization λ',[.00001,.0001,.001,.003,.01,.03,.1,.3,1]),select('method','Algorithm',['ISTA','FISTA'])]
    if family=='compressed':return [slider('P','Measurements P',[256,512,768,1024,1536,2048,3072,4096]),images(('Flower','Cartoon'))]
    raise KeyError(family)


def noise_case(view,p):
    sigma=p['sigma']
    if view in ['noise-example-1d','sparse-signal']:
        t,f=input_signal('Piecewise smooth')
        if view=='sparse-signal':f=np.zeros(256);f[[21,48,106,183,214]]=[1,-.7,.4,-.9,.7]
        y=f+sigma*noise(f.shape)
        return state(plot('Clean signal',series(t,f,'f₀',stems=view=='sparse-signal'),xlabel='t',ylabel='Amplitude'),plot('Noisy observation',series(t,y,'f₀ + w'),series(t,f,'f₀',dash=True),xlabel='t',ylabel='Amplitude'),metrics=[metric('σ',sigma),metric('SNR',snr(f,y),'dB')])
    f=input_image(p['image']);y=f+sigma*noise(f.shape)
    if view=='noise-examples':
        # Variance at unit intensity is sigma² for all three models.
        peak=1/max(sigma*sigma,1e-8);u=np.random.default_rng(SEED).random(f.shape);po=special.pdtr(0,0) if False else None
        rng=np.random.default_rng(SEED);po=rng.poisson(peak*f)/peak if sigma else f;ga=rng.gamma(peak,1/peak,f.shape)*f if sigma else f
        return state(image(y,'Additive Gaussian'),image(po,'Poisson counts / peak'),image(ga,'Multiplicative Gamma'),note='At intensity one the variance equals σ² in each model. Poisson and Gamma noise remain signal-dependent at other intensities.')
    return state(image(f,'Clean f₀'),image(y,'Observation f = f₀ + w'),image(sigma*noise(f.shape),'Noise w',signed=True,scale=.6),metrics=[metric('σ',sigma),metric('SNR',snr(f,y),'dB')],note='w ∼ N(0, σ² I). The same normalized white-noise realization is scaled by σ. Images are clipped only when displayed.')


@functools.lru_cache(maxsize=200)
def gauss_case(name,width):
    f=input_image(name);y=f+.08*noise(f.shape);g=ndimage.gaussian_filter(y,width,mode='wrap') if width else y.copy();return f,y,g

def gaussian(view,p):
    width=p['width'];f,y,g=gauss_case(p['image'],width);t,s=input_signal('Piecewise smooth');sy=s+.08*noise(s.shape);sg=ndimage.gaussian_filter1d(sy,width*2,mode='wrap') if width else sy
    if view in ['iltering-optimal-curve','filtering-optimal-1d']:
        grid=[0,.25,.5,.75,1,1.5,2,3,5];scores=[snr(f,gauss_case(p['image'],v)[2]) for v in grid]
        return state(plot('Filtering risk',series(grid,scores,'SNR'),series([width],[snr(f,g)],'Current width',points=True),xlabel='Gaussian width (pixels)',ylabel='SNR (dB)'),image(g,'Filtered image'),plot('Filtered signal',series(t,s,'Clean'),series(t,sg,'Filtered'),xlabel='t',ylabel='Amplitude'),note='SNR uses the known clean input. The maximum is an oracle choice for this fixed synthetic noise sample.')
    if view=='filtering-denoised':
        n=np.arange(-25,26);kernel=np.exp(-.5*(n/max(2*width,.01))**2);kernel/=kernel.sum()
        return state(plot('Noisy and filtered signal',series(t,sy,'Noisy'),series(t,sg,'Filtered'),xlabel='t',ylabel='Amplitude'),plot('Normalized Gaussian kernel',series(n,kernel,'hσ'),xlabel='Sample offset',ylabel='Weight'),plot('Fourier magnitudes',series(np.arange(129),abs(fft.rfft(sy)),'Noisy'),series(np.arange(129),abs(fft.rfft(sg)),'Filtered'),xlabel='k',ylabel='Magnitude'))
    return state(image(f,'Clean input'),image(y,'Noisy input: σ = 0.08'),image(g,'Gaussian filtered'),metrics=[metric('Width',width,'pixels'),metric('SNR',snr(f,g),'dB')],note='Increasing width reduces noise but blurs edges; the random realization is held fixed.')


@functools.lru_cache(maxsize=1600)
def threshold_case(name,factor,method='Soft',stationary=False,block=1,mu=2):
    f=input_image(name,64);sigma=.08;y=f+sigma*noise(f.shape);fun=C.stationary_denoise if stationary else C.wav_denoise
    mode=method.lower();g=fun(y,factor*sigma,mode,level=3,wavelet='db2',block=block,mu=mu)
    return f,y,g

def threshold(view,p):
    name=p.get('image','Flower');factor=p.get('factor',2);method=p.get('method','Stein' if 'block' in view or view=='comp-thresh' else 'Soft');block=int(p.get('block',1));mu=p.get('mu',2);stationary=view in ['wavthresh-2d-soft-ti','wavthresh-ti','tiwav','tiwav-coefs-thresh','comp-thresh']
    if view=='cycle-spin-principle':
        f=input_image('Flower',64);y=f+.08*noise(f.shape);shift=p['shift'];shifted=np.roll(y,shift,1)
        fun=(lambda z:C.wav_denoise(z,.16,'soft',wavelet='haar',level=3)) if p['method']=='Orthogonal' else (lambda z:C.cycle_denoise(z,.16,'soft',shifts=8,wavelet='haar',level=3))
        a=fun(y);b=np.roll(fun(shifted),-shift,1)
        return state(image(a,'D(f)'),image(b,'S⁻¹ D(Sf)'),image(a-b,'Commutation error',signed=True,scale=.08),metrics=[metric('Translation defect',np.linalg.norm(a-b))],note='The two routes evaluate the shift–denoise–unshift commutative diagram. Averaging all eight Haar phases at three levels restores translation equivariance (up to roundoff).')
    if view in ['semisoft','semisoft-optimal']:method='Semisoft'
    if view=='stein-optimal-curve':method='Stein'
    f,y,g=threshold_case(name,factor,method,stationary,block,mu)
    if view=='semisoft':
        x=np.linspace(-.5,.5,501);T=.08*factor
        return state(plot('Shrinkage rules',series(x,C.semisoft(x,T,mu),'Semi-soft'),series(x,C.stein(x,T),'Stein'),series(x,x,'Identity',dash=True),xlabel='Coefficient a',ylabel='Thresholded value'),image(g,'Semi-soft reconstruction'),metrics=[metric('SNR',snr(f,g),'dB')])
    if view in ['wavthresh-2d-hard-vs-soft','wavthresh-ti','semisoft-optimal','stein-optimal-curve','block-optimal-curves']:
        grid=[.5,1,1.5,2,2.5,3,4];methods=['Hard','Soft'] if view=='wavthresh-2d-hard-vs-soft' else [method];curves=[]
        for mode in methods:curves.append(series(grid,[snr(f,threshold_case(name,r,mode,stationary,block,mu)[2]) for r in grid],mode))
        return state(plot('Threshold selection',*curves,series([factor],[snr(f,g)],'Current selection',points=True),xlabel='T / σ',ylabel='SNR (dB)'),image(g,'Reconstruction'),metrics=[metric('SNR',snr(f,g),'dB')],note='SNR is measured against the known clean image. These are fixed-realization oracle comparisons, not an estimate available from noisy data alone.')
    a,s=C.wavelet_array(y,'db2',3);ag,_=C.wavelet_array(g,'db2',3)
    if view in ['tiwav','tiwav-coefs-thresh']:
        c=pywt.swtn(y,'db2',level=3);z=c[-1]['ad'];q={'Hard':C.hard,'Soft':C.soft,'Stein':C.stein}[method](z,.08*factor)
        return state(image(y,'Noisy image'),image(z,'Undecimated horizontal details',signed=True,scale=.4),image(q,'Thresholded details',signed=True,scale=.4),note='Undecimated coefficients are full-size images at every scale. Their signs share a fixed display scale.')
    return state(image(y,'Noisy input'),waveimage(ag,s,'Reconstructed wavelet coefficients'),image(g,method+' reconstruction'),metrics=[metric('SNR',snr(f,g),'dB'),metric('Threshold',factor*.08),metric('Block side',block)],note=('Translation-invariant' if stationary else 'Orthogonal')+' db2 wavelets with three levels; coarse coefficients are retained. Block Stein shrinkage uses the mean squared coefficient magnitude in each block.')


def poisson(view,p):
    from scipy.stats import poisson as distribution
    lam=p['photons'];f=input_image(p['image'],64);rng=np.random.default_rng(SEED);u=rng.random(f.shape);counts=distribution.ppf(u,lam*f);y=counts/lam
    if view=='poisson-distrib':
        k=np.arange(max(15,int(lam+5*np.sqrt(lam))));probs=distribution.pmf(k,lam)
        return state(plot('Poisson distribution',series(k,probs,'P(X = k)',bars=True),xlabel='Count k',ylabel='Probability'),metrics=[metric('Mean',lam),metric('Variance',lam)])
    if view=='variance-stabilization':
        grid=np.geomspace(.1,128,60);curves=[]
        for name,fn in [('Anscombe',lambda k:2*np.sqrt(k+3/8)),('Freeman–Tukey',lambda k:np.sqrt(k)+np.sqrt(k+1)),('Square root',lambda k:2*np.sqrt(k))]:
            variances=[]
            for v in grid:
                k=np.arange(max(30,int(v+12*np.sqrt(v)+30)));prob=distribution.pmf(k,v);z=fn(k);variances.append(np.dot(prob,z*z)-np.dot(prob,z)**2)
            curves.append(series(grid,variances,name))
        return state(plot('Exact transformed variances',*curves,xlabel='Poisson mean',ylabel='Variance',logx=True),image(y,'Poisson observation'),note='Variances are summed from the Poisson probability mass function; stabilization approaches variance one at sufficiently large means.')
    if view=='poisson-wav':
        raw=C.wav_denoise(y,2/np.sqrt(lam),'soft',wavelet='db2',level=3);z=2*np.sqrt(counts+3/8);d=C.wav_denoise(z,2,'soft',wavelet='db2',level=3);g=np.maximum((d/2)**2-3/8,0)/lam
        return state(image(y,'Poisson observation'),image(raw,'Without stabilization'),image(g,'Anscombe stabilization'),metrics=[metric('Direct SNR',snr(f,raw),'dB'),metric('Stabilized SNR',snr(f,g),'dB')],note='The algebraic Anscombe inverse is used here; it is biased at very low counts. The direct method uses the worst-case standard deviation 1/√λmax.')
    return state(image(f,'Clean intensity f₀'),image(y,'Counts / λmax'),metrics=[metric('Peak photons',lam),metric('SNR',snr(f,y),'dB')],note='N(x) ∼ Poisson(λmax f₀(x)). A fixed uniform sample is transformed through each Poisson inverse CDF, keeping slider comparisons coupled.')


def multiplicative(view,p):
    from scipy.stats import gamma
    K=p['looks'];f=input_image(p['image'],64);u=np.random.default_rng(SEED).random(f.shape);eta=gamma.ppf(u,K,scale=1/K);y=f*eta;bias=special.digamma(K)-np.log(K);sd=np.sqrt(special.polygamma(1,K));z=np.log(np.maximum(y,1e-8))-bias
    if view=='multnoise-hist-unstab':
        raw,edges=np.histogram(eta,bins=40,density=True);st,ed=np.histogram((np.log(eta)-bias)/sd,bins=40,density=True)
        return state(plot('Multiplicative noise',series((edges[1:]+edges[:-1])/2,raw,'Gamma density',bars=True),xlabel='η',ylabel='Density'),plot('Centered standardized log noise',series((ed[1:]+ed[:-1])/2,st,'Log noise',bars=True),xlabel='(log η − bias) / sd',ylabel='Density'),metrics=[metric('Variance of η',1/K)])
    if view=='multnoise-wav':
        d=C.wav_denoise(z,2*sd,'soft',wavelet='db2',level=3);g=np.exp(d);raw=C.wav_denoise(y,2/np.sqrt(K),'soft',wavelet='db2',level=3)
        return state(image(y,'Gamma-corrupted image'),image(g,'Log-domain thresholding'),image(raw,'Direct thresholding'),metrics=[metric('Log-domain SNR',snr(f,g),'dB'),metric('Direct SNR',snr(f,raw),'dB')],note='The log-noise mean ψ(K) − log K is removed and its exact variance ψ₁(K) sets the threshold. Exponentiation is not an unbiased intensity estimator.')
    return state(image(f,'Clean image'),image(y,'Multiplicative Gamma observation'),metrics=[metric('Looks K',K),metric('Relative standard deviation',1/np.sqrt(K))],note='η ∼ Gamma(K, scale = 1/K), so Eη = 1 and Varη = 1/K.')


def gradient(view,p):
    f=input_image(p['image']);eps=p.get('eps',.01);dy,dx=C.gradient(f);norm=np.sqrt(dx*dx+dy*dy+eps*eps);scale=np.quantile(abs(np.r_[dx.ravel(),dy.ravel()]),.98);gain=1/p.get('gain',1)
    def tint(z,color):
        # Keep the same hue on either side of zero. Varying just one RGB
        # channel would turn negative red into cyan and negative blue into yellow.
        a=np.sign(z)*np.power(np.minimum(abs(z)/(scale*gain+1e-9),1),.65)
        color=np.array(color);positive=.05+.95*color;negative=.15*color
        target=np.where(a[...,None]>=0,positive,negative)
        return .5+abs(a[...,None])*(target-.5)
    if view=='grad-discrete':
        coarse=f[::8,::8];gy,gx=C.gradient(coarse);panel=image(coarse,'Forward differences on the pixel grid');boost=p.get('gain',1)
        panel['vectors']=[{'x':float(x)+.5,'y':float(y)+.5,'dx':float(gx[y,x])*boost*3,'dy':float(gy[y,x])*boost*3} for y in range(16) for x in range(16)]
        return state(panel,image(np.sqrt(gx*gx+gy*gy),'Gradient magnitude'),metrics=[metric('Arrow display gain',boost)],note='Arrows show (dx, dy) from periodic forward differences on a 16 × 16 grid. The vertical coordinate increases downward, as in the image. Gain scales arrows only.')
    if view=='discrete-operators-lapl':return state(image(C.laplacian(f),'Discrete Laplacian Δf',signed=True),image(C.divergence(C.gradient(f)/norm),'−∇TVε(f)',signed=True),note='Periodic forward differences and their negative adjoint divergence are used. The TV force divides the gradient by √(|∇f|² + ε²).')
    if view=='discrete-operators-tveps':return state(image(f,'Input'),image(norm,'√(|∇f|² + ε²)'),metrics=[metric('Regularized TV',norm.sum()),metric('ε',eps)])
    red=image(tint(dx,[1,0,0]),'Horizontal derivative dx — red');blue=image(tint(dy,[0,0,1]),'Vertical derivative dy — blue')
    for panel,negative,positive in [(red,'#260000','#ff0d0d'),(blue,'#000026','#0d0dff')]:
        panel['colorLabels']=[{'label':'Negative','color':negative},{'label':'Zero','color':'#808080'},{'label':'Positive','color':positive}]
    return state(image(f,'Input f'),red,blue,metrics=[metric('Display gain',p.get('gain',1))],note='Neutral gray means zero; brighter red or blue is positive, darker is negative. A shared 98th-percentile scale and a 0.65 power enhance visibility. The gain changes display contrast, not the finite differences.')


@functools.lru_cache(maxsize=400)
def pde_case(name,method,strength,noisy):
    f=input_image(name,64);y=f+.08*noise(f.shape) if noisy else f
    if not strength:return f,y,y.copy()
    if method=='Heat flow':g=C.heat(y,8*strength)
    elif method=='TV flow':g=tv_flow(y,[strength],eps=.02)[0]
    elif method=='Sobolev regularization':g=C.sobolev_denoise(y,8*strength)
    else:g=C.tv_denoise(y,strength,300)
    return f,y,g

def pde(view,p):
    strength=p['strength'];method=p['method'];noisy=view!='flow';f,y,g=pde_case(p['image'],method,strength,noisy)
    if view=='denoising-cameraman-snr':
        vals=[0,.01,.03,.06,.1,.2,.4];curves=[]
        for m in ['Heat flow','TV flow','Sobolev regularization','TV regularization']:curves.append(series(vals,[snr(f,pde_case(p['image'],m,v,True)[2]) for v in vals],m))
        return state(plot('Denoising parameter sweep',*curves,xlabel='Control strength',ylabel='SNR (dB)'),image(g,method),note='For comparable visible ranges, heat time and Sobolev λ are eight times the control strength; TV time and TV λ equal it. Curves are oracle SNR values.')
    return state(image(y,'Noisy input' if noisy else 'Input'),image(g,method),image(g-y,'Change from input',signed=True),metrics=[metric('t or λ',strength*(8 if method in ['Heat flow','Sobolev regularization'] else 1)),metric('SNR',snr(f,g),'dB'),metric('Discrete TV',C.tv(g))],note='Periodic finite differences. Heat flow uses its exact Fourier multiplier; TV flow uses stable explicit steps; regularization solves a fidelity-plus-prior objective.')


@functools.lru_cache(maxsize=100)
def inpaint_case(steps,method,maskname):
    f=input_image('Flower',64);n=64;y,x=np.mgrid[:n,:n]
    mask=~((x>20)&(x<44)&(y>20)&(y<44)) if maskname=='Square' else ((x%12>3)&(y%12>3) if maskname=='Stripes' else np.random.default_rng(SEED).random(f.shape)>.6)
    observed=f*mask
    if steps==0:g=observed.copy()
    elif method in ['Sobolev','TV']:g=C.inpaint(observed,mask,kind=method.lower(),steps=steps)
    else:
        W,S=wavelet_operator(f.shape);A=lambda a:mask*S(a);At=lambda v:W(mask*v);g=S(ista(observed,A,At,.01,steps=steps,accelerated=True))
    return f,mask,observed,g

def inpainting(view,p):
    f,mask,y,g=inpaint_case(p['steps'],p['method'],p['mask']);display=np.where(mask,y,1)
    return state(image(f,'Reference'),image(display,'Observed pixels (white = missing)'),image(g,p['method']+' reconstruction'),metrics=[metric('Observed fraction',mask.mean()),metric('SNR',snr(f,g),'dB'),metric('Iterations',p['steps'])],note='Sobolev and TV inpainting enforce the known samples exactly. Wavelet inpainting minimizes a soft fidelity term plus 0.01‖a‖₁ using FISTA; all methods use the same mask.')


@functools.lru_cache(maxsize=30)
def tomography_case(lines,method):
    f=C.phantom(64);mask=radial_mask(64,lines);y=fft.fft2(f,norm='ortho')*mask;g=fourier_tv(y,mask,.02,600) if method=='TV' else fft.ifft2(y,norm='ortho').real
    return f,mask,g,y

def tomography(view,p):
    if view=='tomo-principle':
        theta=np.deg2rad(p['angle']);t=np.linspace(-1.2,1.2,401);proj=C.phantom_projection(t,theta)
        return state(image(C.phantom(128),'Modified Shepp–Logan phantom'),plot('Exact Radon projection',series(t,proj,'Rf(θ, t)'),xlabel='Detector coordinate t',ylabel='Line integral'),metrics=[metric('θ',p['angle'],'°')],note='The projection integrates the same ten analytic ellipses used to construct the phantom; θ is the detector-normal angle.')
    f,mask,g,y=tomography_case(p['lines'],p['method']);return state(image(f,'Phantom'),image(fft.fftshift(mask).astype(float),'Measured Fourier samples'),image(g,p['method']+' reconstruction'),metrics=[metric('Sample fraction',mask.mean()),metric('SNR',snr(f,g),'dB'),metric('Relative data residual',np.linalg.norm(mask*fft.fft2(g,norm='ortho')-y)/np.linalg.norm(y))],note='This is the book’s Cartesian radial Fourier sampling model, not filtered backprojection from a sinogram. The mask is Hermitian; TV reconstruction uses λ = 0.02 and 600 primal–dual steps.')


@functools.lru_cache(maxsize=100)
def sparse_case(lam,method,oned=False):
    if oned:
        n=128;x=np.zeros(n);x[[14,33,61,83,109]]=[1,-.8,.7,1,-.6]
        if oned=='signal':x=input_signal('Piecewise smooth',n)[1]
        H=blur_symbol(x.shape,1.3);y=convolution(x,H)+.01*noise(x.shape)
        if oned=='signal':
            coeff=pywt.wavedec(x,'db2',level=3,mode='periodization');cuts=np.cumsum([len(v) for v in coeff])[:-1]
            W=lambda f:np.concatenate(pywt.wavedec(f,'db2',level=3,mode='periodization'));S=lambda a:pywt.waverec(np.split(a,cuts),'db2',mode='periodization')
        else:W=lambda x:x;S=lambda x:x
        A=lambda c:convolution(S(c),H);At=lambda z:W(convolution(z,H));a,E,it=ista(y,A,At,lam,steps=240,accelerated=method=='FISTA',history=True);return x,y,S(a),E,it
    f=input_image('Flower',64);H=blur_symbol(f.shape,1.2);y=convolution(f,H)+.015*noise(f.shape);W,S=wavelet_operator(f.shape);A=lambda c:convolution(S(c),H);At=lambda z:W(convolution(z,H));c,E,it=ista(y,A,At,lam,steps=240,accelerated=method=='FISTA',history=True);return f,y,S(c),E,None

def sparse(view,p):
    if view=='sparsity-wav':
        f=input_image(p['image']);a,s=C.wavelet_array(f,'db2',3);T=p['factor'];b=C.hard(a,T);ordered=np.sort(abs(a).ravel())[::-1]
        return state(image(f,'Image'),waveimage(b,s,'Coefficients above threshold'),plot('Ordered magnitudes',series(np.arange(1,len(ordered)+1)[::20],ordered[::20],'|a|(m)'),xlabel='Rank m',ylabel='Magnitude',logx=True),metrics=[metric('Nonzeros',np.count_nonzero(b))])
    lam=p['lam'];method=p['method'];one=view in ['spikes-lun','spikes-ista-decay','deconvolution-1d'];f,y,g,E,it=sparse_case(lam,method,'signal' if view=='deconvolution-1d' else one)
    if view in ['spikes-ista-decay','deconvolution-2d-snr']:
        if view=='spikes-ista-decay':return state(plot('Objective along the iterates',series(np.arange(1,len(E)+1),E,method),xlabel='Iteration',ylabel='Objective',logy=True),plot('Sparse estimate',series(np.arange(len(f)),f,'Truth',stems=True),series(np.arange(len(f)),g,'Estimate',stems=True),xlabel='Index',ylabel='Coefficient'))
        vals=[.00001,.0001,.001,.003,.01,.03,.1,.3,1];scores=[snr(f,sparse_case(v,method)[2]) for v in vals]
        return state(plot('Regularization parameter',series(vals,scores,method),series([lam],[snr(f,g)],'Current λ',points=True),xlabel='λ',ylabel='SNR (dB)',logx=True),image(g,'Wavelet-sparse deconvolution'))
    if one:return state(plot('Observation',series(np.arange(len(f)),f,'Truth',stems=True),series(np.arange(len(y)),y,'Blurred and noisy'),xlabel='Sample',ylabel='Value'),plot('Sparse deconvolution',series(np.arange(len(f)),f,'Truth',stems=True),series(np.arange(len(g)),g,'Estimate',stems=True),xlabel='Sample',ylabel='Coefficient'),metrics=[metric('λ',lam),metric('SNR',snr(f,g),'dB')],note=('A piecewise-smooth signal is reconstructed in an orthogonal db2 wavelet basis.' if view=='deconvolution-1d' else 'Sparse spikes are reconstructed in the coordinate basis.')+' A norm-one periodic Gaussian blur is inverted with ℓ¹ regularization.')
    return state(image(f,'Reference'),image(y,'Blurred and noisy'),image(g,'Wavelet-sparse deconvolution'),metrics=[metric('λ',lam),metric('SNR',snr(f,g),'dB'),metric('Final objective',E[-1])],note='Minimizes ½‖A Wa − y‖² + λ‖a‖₁ with an orthogonal db2 synthesis operator and a norm-one blur; 240 iterations.')


def hadamard(v):
    z=np.array(v,copy=True);h=1
    while h<len(z):
        b=z.reshape(-1,2*h);left=b[:,:h].copy();right=b[:,h:].copy();b[:,:h]=left+right;b[:,h:]=left-right;h*=2
    return z/np.sqrt(len(z))

@functools.lru_cache(maxsize=30)
def compressed_case(name,P):
    f=input_image(name,64);N=f.size;rng=np.random.default_rng(SEED);order=np.r_[0,rng.permutation(np.arange(1,N))];scramble=rng.permutation(N);rows=order[:P]
    def A(x):return hadamard(x.ravel()[scramble])[rows]
    def At(y):
        z=np.zeros(N);z[rows]=y;out=np.empty(N);out[scramble]=hadamard(z);return out.reshape(f.shape)
    W,S=wavelet_operator(f.shape);y=A(f);c=ista(y,lambda c:A(S(c)),lambda y:W(At(y)),.0005,steps=1200,accelerated=True);g=S(c)
    probe=rng.standard_normal(f.shape);v=rng.standard_normal(P);assert abs(np.dot(A(probe),v)-np.sum(probe*At(v)))<1e-9
    return f,At(y),g,np.linalg.norm(A(g)-y)/np.linalg.norm(y)

def compressed(view,p):
    f,back,g,res=compressed_case(p['image'],p['P']);return state(image(f,'Reference: 4096 pixels'),image(back,'Adjoint backprojection'),image(g,'Wavelet ℓ¹ reconstruction'),metrics=[metric('Measurements P',p['P']),metric('P / Q',p['P']/4096),metric('SNR',snr(f,g),'dB'),metric('Relative data residual',res)],note='Simulated signed Walsh–Hadamard measurements with the same fixed pixel permutation and nested row sets. FISTA minimizes the actual measurement fidelity plus 0.0005‖a‖₁ for 1200 iterations.')

FUNCTIONS={'noise':noise_case,**{k:globals()[k] for k in ['gaussian','threshold','poisson','multiplicative','gradient','pde','inpainting','tomography','sparse','compressed']}}
