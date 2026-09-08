"""Sampling, representation, source coding, and approximation experiments."""
from numerics import *
import common as C
from scipy.spatial import Delaunay

def controls(family,view):
    if family=='acquisition':
        if view=='examples-2':return [slider('channel','Spectral channel',[1,4,8,12,16,20,24,28,32])]
        return [slider('zoom','Zoom factor',[1,2,4,8,16]),select('source','Source',['Flower','Birdsong'])]
    if family=='quantization':return [slider('K','Gray levels K',[2,3,4,8,16,32,64]),images()]
    if family=='fourier':
        if view=='filtering-1d':return [slider('sigma','Filter width σ (samples)',[0,.5,1,2,4,8,16]),signals()]
        if view=='2d-fourier-extension':return [slider('k','Horizontal frequency',[0,1,2,4,6,8]),slider('l','Vertical frequency',[0,1,2,4,6,8])]
        if view=='fourier-2d':return [slider('taper','Boundary taper',[0,.2,.4,.6,.8,1]),images()]
        return [slider('scale','Scale level',[1,2,3,4]),select('wavelet','Wavelet',['haar','db2','db4'])]
    if family=='coding':return [slider('skew','Distribution concentration',[0,.25,.5,1,2,3,4]),slider('K','Alphabet size',[2,3,4,6,8])]
    if family=='wavelet1d':return [slider('level','Decomposition level',[1,2,3,4]),select('wavelet','Wavelet',['haar','db2','db4']),signals()]
    if family=='wavelet2d':return [slider('level','Decomposition level',[1,2,3,4]),select('wavelet','Wavelet',['haar','db2','db4']),images(('Flower','Cartoon'))]
    if family=='models':
        if view=='image-model-cartoon':return [slider('count','Opaque leaves',[5,15,40,100,250,600])]
        return [slider('strength','Smoothing strength',[0,.5,1,2,4,8]),images(('Flower','Cartoon','Mandrill','Film cartoon'))]
    if family=='approximation1d':return [slider('M','Retained coefficients M',[4,8,16,32,64,128,256]),select('basis','Basis',['Fourier','haar','db4']),signals()]
    if family=='approximation2d':return [slider('M','Retained coefficients M',[16,64,256,1024,4096,8192,16384]),select('basis','Basis',['Fourier','DCT','haar','db4']),images()]
    if family=='mesh':return [slider('angular','Vertices around the boundary',[32,40,48,56,64]),slider('variation','Smooth variation',[.4,.7,1,1.3])]
    if family=='directional':
        if view=='curvelet-denoising':return [slider('factor','Threshold T / σ',[1,1.5,2,2.5,3,3.5,4]),images(('Flower','Cartoon'))]
        return [slider('scale','Radial scale',[4,8,16,32]),slider('angle','Orientation (degrees)',[0,30,60,90,120,150])]
    if family=='compression':return [slider('rate','Target bits / pixel',[.3,.5,.75,1,1.5,2,3]),images(('Flower','Mandrill'))]
    raise KeyError(family)


def acquisition(view,p):
    if view=='examples-2':
        f=C.flower(128,color=True);k=int(p['channel'])-1;w=400+k*300/31
        centers=np.array([620,540,460]);weights=np.exp(-.5*((w-centers)/40)**2);z=f@weights
        spec=np.array([f[64,64]@np.exp(-.5*((v-centers)/40)**2) for v in np.linspace(400,700,32)])
        return state(image(f,'RGB flower'),image(z,f'Synthetic channel: {w:.0f} nm'),plot('Spectrum at the central pixel',series(np.linspace(400,700,32),spec,'Intensity'),series([w],[spec[k]],'Selected channel',points=True),xlabel='Wavelength (nm)',ylabel='Intensity'),note='The 32 channels are a synthetic spectral model derived from RGB, not measured multispectral data.')
    z=p['zoom']
    if p['source']=='Birdsong':
        fs,y=C.bird();t=np.arange(len(y))/fs;a=int(len(y)*.4);count=int(len(y)/z*.4);ids=np.linspace(0,len(y)-1,1000).astype(int);ids2=np.linspace(a,min(a+count,len(y)-1),900).astype(int)
        return state(plot('Recorded birdsong',series(t[ids],y[ids],'Amplitude'),xlabel='Time (s)',ylabel='Amplitude'),plot('Short-time detail',series(t[ids2],y[ids2],'Amplitude'),xlabel='Time (s)',ylabel='Amplitude'),metrics=[metric('Sampling rate',fs,'Hz')],note='The waveform comes from data/bird.wav; the audio player contains the same recording.')
    f=C.flower(128,color=True);size=128//z;start=(128-size)//2;box={'x':start,'y':start,'w':size,'h':size,'color':'#cc3333'}
    return state(image(f,'Flower — selected region',boxes=[box]),image(f[start:start+size,start:start+size],'Pixelated enlargement',boxes=[{'x':0,'y':0,'w':size,'h':size,'color':'#cc3333'}]),metrics=[metric('Crop width',size,'pixels')])


def quantization(view,p):
    f=input_image(p['image']);K=p['K'];g=(np.minimum(np.floor(K*f),K-1)+.5)/K
    x=np.linspace(0,1,501);q=(np.minimum(np.floor(K*x),K-1)+.5)/K
    return state(image(f,'Input f'),image(g,f'Midpoint quantization: K = {K}'),plot('Quantizer',series(x,q,'QK(x)'),series(x,x,'Identity',dash=True),xlabel='Input intensity',ylabel='Reconstruction',ylim=[0,1]),metrics=[metric('SNR',snr(f,g),'dB'),metric('Maximum error bound',.5/K)],note='Uniform bins on [0,1], with endpoint clipping. Values are reconstructed at each bin midpoint.')


def fourier(view,p):
    if view=='filtering-1d':
        t,f=input_signal(p['signal']);s=p['sigma'];g=ndimage.gaussian_filter1d(f,s,mode='wrap') if s else f.copy()
        return state(plot('Gaussian filtering',series(t,f,'Input'),series(t,g,'Filtered'),xlabel='t',ylabel='f(t)'),plot('Fourier magnitudes',series(np.arange(129),abs(fft.rfft(f)),'Input'),series(np.arange(129),abs(fft.rfft(g)),'Filtered'),xlabel='Frequency k',ylabel='Magnitude'),metrics=[metric('Filter width',s,'samples')])
    if view=='2d-fourier-extension':
        y,x=np.mgrid[:128,:128]/128;k=p['k'];l=p['l'];z=np.exp(2j*np.pi*(k*x+l*y))
        return state(image(z.real,'Real part: cos(2π(kx + ℓy))',signed=True),image(z.imag,'Imaginary part: sin(2π(kx + ℓy))',signed=True),metrics=[metric('Horizontal frequency',k),metric('Vertical frequency',l)])
    if view=='fourier-2d':
        f=input_image(p['image']);w=signal.windows.tukey(128,alpha=p['taper']);g=f*w[:,None]*w[None,:]
        def spectrum(z):return np.log1p(abs(fft.fftshift(fft.fft2(z))))/np.log1p(np.sum(f))
        return state(image(f,'Original'),image(g,'Windowed image'),image(spectrum(f),'Original log spectrum'),image(spectrum(g),'Windowed log spectrum'),note='A separable Tukey window suppresses the boundary jump introduced by periodic extension. The two spectra share one display scale.')
    t=np.linspace(-.1,1.1,700);L=p['scale'];phi,psi,x=pywt.Wavelet(p['wavelet']).wavefun(level=8);t=np.linspace(-.1,(x[-1]+3)/2**L,700)
    curves=[series(t,np.interp(2**L*t-k,x,psi,left=0,right=0)*2**(L/2),f'n = {k}') for k in range(min(2**L,4))]
    return state(plot('Fourier modes',*[series(t,np.sqrt(2)*np.cos(2*np.pi*k*t),f'k = {k}') for k in [1,2,4]],xlabel='t',ylabel='Atom value'),plot('Translated wavelets',*curves,xlabel='t',ylabel='Scaled atom'),note='Translations share the selected dyadic scale; Fourier atoms remain global. The wavelet uses its canonical support and the normalization 2^(L/2).')


def coding(view,p):
    K=p['K'];skew=p.get('skew',1);prob=np.exp(-skew*np.arange(K));prob/=prob.sum();codes=C.huffman_codes(dict(enumerate(prob)))
    H=-np.sum(prob*np.log2(prob));L=sum(prob[k]*len(c) for k,c in codes.items());rows=[[f's{k+1}',f'{prob[k]:.4f}',codes[k],len(codes[k])] for k in range(K)]
    assert H<=L+1e-10 and L<H+1+1e-10
    dist=plot('Symbol probabilities',series(np.arange(1,K+1),prob,'pk',bars=True),xlabel='Symbol k',ylabel='Probability',ylim=[0,1])
    if view=='the-concave-function-underlying-entropy':
        x=np.linspace(.0001,1,401);return state(plot('Entropy contribution h(u)',series(x,-x*np.log2(x),'−u log₂ u'),series(prob,-prob*np.log2(prob),'Current probabilities',points=True),xlabel='u',ylabel='h(u) (bits)'),dist,metrics=[metric('Entropy H(p)',H,'bits')])
    if view=='codage-differences':
        t,f=input_signal('Smooth');q=np.round(f*(K-1)).astype(int);d=np.diff(q);u,c=np.unique(d,return_counts=True);r=c/c.sum();Hd=-sum(r*np.log2(r));uq,cq=np.unique(q,return_counts=True);rq=cq/cq.sum()
        return state(plot('Quantized signal and differences',series(t,q,'Quantized signal'),series(t[1:],d,'First differences'),xlabel='t',ylabel='Symbol value'),plot('Empirical distributions',series(uq,rq,'Signal',bars=True),series(u,r,'Differences',points=True),xlabel='Symbol',ylabel='Probability'),metrics=[metric('Difference entropy',Hd,'bits'),metric('Original entropy',-sum(rq*np.log2(rq)),'bits')],note='Differencing is invertible when the initial sample is retained. Histograms ignore remaining temporal dependencies.')
    if view=='huff-code':
        clusters=[([str(k+1)],float(v)) for k,v in enumerate(prob)];merges=[]
        while len(clusters)>1:
            clusters.sort(key=lambda v:(v[1],v[0]));a,b=clusters[:2];c=(a[0]+b[0],a[1]+b[1]);merges.append([' + '.join(a[0]),' + '.join(b[0]),f'{c[1]:.4f}']);clusters=[*clusters[2:],c]
        return state(table('Huffman merges',['Cluster 1','Cluster 2','Merged weight'],merges),table('Resulting prefix code',['Symbol','Probability','Word','Length'],rows),metrics=[metric('Entropy',H,'bits'),metric('Expected length',L,'bits')])
    return state(dist,table('Huffman prefix code',['Symbol','Probability','Word','Length'],rows),metrics=[metric('Entropy H(p)',H,'bits'),metric('Expected length L',L,'bits'),metric('Kraft sum',sum(2**-len(c) for c in codes.values()))],note='Words and lengths are recomputed by Huffman’s algorithm. The code has no leaf that is a prefix of another.')


def wavelet1d(view,p):
    t,f=input_signal(p['signal']);w=p['wavelet'];level=p['level'];c,a=wave1(f,w,level);g=wav1_approx(f,w,level);recon=pywt.waverec(c,w,mode='periodization');err=np.linalg.norm(recon-f)
    assert err<1e-10 and abs(sum(np.sum(v*v) for v in c)-np.sum(f*f))<1e-9
    base=plot('Input and approximation',series(t,f,'f'),series(t,g,'PV f'),xlabel='t',ylabel='Amplitude')
    coeff=plot('Packed coefficients',series(np.arange(len(a)),a,'a, d',stems=True),xlabel='Packed index',ylabel='Coefficient')
    if view in ['haard-scaling-wav','wavelets']:
        phi,psi,x=pywt.Wavelet(w).wavefun(level=8);scale=2.**(-level);return state(plot('Scaling function',series(x*scale,phi/np.sqrt(scale),'φj'),xlabel='t',ylabel='Amplitude'),plot('Mother wavelet',series(x*scale,psi/np.sqrt(scale),'ψj'),xlabel='t',ylabel='Amplitude'),metrics=[metric('Vanishing moments',pywt.Wavelet(w).vanishing_moments_psi)],note='Here j = −level, consistent with the book’s convention: larger j means coarser resolution.')
    if view in ['inversion-transpose','haar-analysis-synthesis-matrix']:
        n=32;lev=min(level,2);W=np.stack([wave1(e,w,lev)[1] for e in np.eye(n)],axis=1);assert np.linalg.norm(W.T@W-np.eye(n))<1e-10
        return state(image(W,'Analysis W',signed=True),image(W.T,'Synthesis Wᵀ',signed=True),image(W.T@W,'Wᵀ W = I'),metrics=[metric('Orthogonality error',np.linalg.norm(W.T@W-np.eye(n)))],note='Matrices act on 32 samples, with periodic boundary extension. Columns are transforms of coordinate vectors.')
    if view in ['haar-appr-details','embeded-spaces-1d']:
        coarse=wav1_approx(f,w,level+1);detail=g-coarse
        return state(base,plot('Orthogonal detail',series(t,detail,'PVj f − PVj+1 f'),xlabel='t',ylabel='Detail'),metrics=[metric('Energy in detail',np.sum(detail**2))],note='The finer approximation equals the next coarser approximation plus its orthogonal detail.')
    if view in ['filterbank-fwd','filterbank-pyramidal','filterbank-bwd']:
        return state(base,plot('Approximation vector',series(np.arange(len(c[0])),c[0],'aj',stems=True),xlabel='n',ylabel='aj[n]'),plot('Detail vector',series(np.arange(len(c[1])),c[1],'dj',stems=True),xlabel='n',ylabel='dj[n]'),metrics=[metric('Reconstruction error',err)],note='Filtering and decimation are evaluated with the selected orthogonal wavelet; synthesis reconstructs the original samples.')
    return state(base,coeff,metrics=[metric('Reconstruction error',err),metric('Coarse coefficients',len(c[0]))],note='Periodized orthogonal transform. Packing is coarse approximation first, then details from coarse to fine.')


def wavelet2d(view,p):
    f=input_image(p['image']);w=p['wavelet'];L=p['level'];a,s=C.wavelet_array(f,w,L);recon=C.wavelet_inverse(a,s,w);assert np.linalg.norm(recon-f)<1e-9
    coarse=np.zeros_like(a);coarse[s[0]]=a[s[0]];g=C.wavelet_inverse(coarse,s,w)
    if view=='wavelets-2d':
        phi,psi,x=pywt.Wavelet(w).wavefun(level=5)
        # Downsample only the rendering; the tensor construction is exact.
        i=np.linspace(0,len(phi)-1,128).astype(int);ph=phi[i];ps=psi[i]
        return state(image(np.outer(ph,ps),'Horizontal variation φ ⊗ ψ',signed=True),image(np.outer(ps,ph),'Vertical variation ψ ⊗ φ',signed=True),image(np.outer(ps,ps),'Diagonal variation ψ ⊗ ψ',signed=True),note='Tensor products of the selected one-dimensional scaling and wavelet functions.')
    if view in ['wavalgo-2d-aniso','wavecoefs-iso']:
        row=np.apply_along_axis(lambda v:wave1(v,w,L)[1],1,f);tensor=np.apply_along_axis(lambda v:wave1(v,w,L)[1],0,row);bounds=[]
        for l in range(1,L+1):
            color=C.SCALE_COLORS[(l-1)%len(C.SCALE_COLORS)];bounds.extend([{'x':0,'y':128/2**l,'w':128,'h':0,'color':color},{'x':128/2**l,'y':0,'w':0,'h':128,'color':color}])
        z=.5+.5*np.sign(tensor)*np.minimum(np.log1p(20*abs(tensor))/4,1);z[:128//2**L,:128//2**L]=tensor[:128//2**L,:128//2**L]/2**L
        return state(image(f,'Input'),image(z,'Full tensor transform',boundaries=bounds),waveimage(a,s,'Isotropic pyramid'),note='The tensor transform decomposes every row and column independently. The isotropic pyramid recursively decomposes only the coarse quadrant.')
    if view in ['multires-2d-haar','embeded-spaces-2d']:
        return state(image(f,'Input f'),image(g,'Projection PV f'),image(f-g,'Omitted details',signed=True),metrics=[metric('Retained coarse coefficients',coarse[s[0]].size),metric('Approximation SNR',snr(f,g),'dB')])
    if view in ['wavalgo-2d-step-fwd','wavalgo-2d-step-bwd','wavalgo-2d']:
        c=pywt.dwt2(f,w,mode='periodization');parts=[image(c[0]/2,'a: coarse average')]+[image(v,'d: '+name,signed=True) for v,name in zip(c[1],['vertical','horizontal','diagonal'])]
        return state(*parts,metrics=[metric('One-step inverse error',np.linalg.norm(pywt.idwt2(c,w,mode='periodization')-f))],note='One separable analysis step: two low-pass filters give a; at least one high-pass filter gives each detail subband.')
    return state(image(f,'Input image'),waveimage(a,s),image(g,'Coarse reconstruction'),metrics=[metric('Parseval residual',abs(np.sum(a*a)-np.sum(f*f)))],note='Colored lines mark actual dyadic subband boundaries. Signed logarithmic contrast affects only the coefficient display.')


def models(view,p):
    if view=='image-model-cartoon':
        rng=np.random.default_rng(SEED);y,x=np.mgrid[:128,:128]/128;f=np.ones((128,128))*.85
        for _ in range(p['count']):
            cx,cy=rng.random(2);r=rng.uniform(.015,.22);f[(x-cx)**2+(y-cy)**2<r*r]=rng.uniform(.08,.92)
        return state(image(f,'Dead-leaves construction'),metrics=[metric('Opaque leaves',p['count'])],note='The same seeded sequence of opaque discs is used throughout. Later discs occlude earlier ones.')
    strength=p['strength'];white=noise((128,128));f=input_image(p['image']);s=ndimage.gaussian_filter(white,strength,mode='wrap') if strength else white
    s=(s-s.mean())/s.std()*.15+.5
    if view=='approximation-class':
        tv=C.tv_denoise(.5+.15*white,strength*.05,250) if strength else .5+.15*white
        return state(image(input_image('Smooth'),'Smooth model'),image(tv,'BV: TV-denoised white noise'),image(input_image('Cartoon'),'Cartoon model'),note='The bounded-variation example solves the discrete ROF problem on a pure white-noise image. Display clipping does not change the computation.')
    if view=='signal-model-sobolev':return state(image(.5+.15*white,'White noise'),image(s,'Gaussian-smoothed noise'),note='Mean and variance are normalized equally for display, so changes in correlation length remain visible.')
    if view=='approx-speed-1':return state(image(f,p['image']),image(ndimage.gaussian_filter(f,strength,mode='wrap') if strength else f,'Smoothed image'),metrics=[metric('Smoothing width',strength,'pixels')])
    if view in ['sample-cartoon','singular-part']:
        g=ndimage.gaussian_filter(f,strength,mode='wrap') if strength else f
        return state(image(f,p['image']),image(g,'Smoothed model'),image(np.sqrt(np.sum(C.gradient(g)**2,axis=0)),'Gradient magnitude'),note='The film image, when selected, is the documented public-domain Feline Follies frame (1919).')
    raise KeyError(view)


def approximation1d(view,p):
    from approximation_models import real_fourier
    t,f=input_signal(p['signal']);M=p['M'];w=p['basis']
    if w=='Fourier':a,inv,r=real_fourier(f);order=np.argsort(r,kind='stable')
    else:
        c,a=wave1(f,w,4);lengths=[len(v) for v in c];cuts=np.cumsum(lengths)[:-1];inv=lambda b:pywt.waverec(np.split(b,cuts),w,mode='periodization');order=np.arange(len(a))
    b=np.zeros_like(a);b[order[:M]]=a[order[:M]];linear=inv(b);s=C.keep_largest(a,M);nonlin=inv(s);assert np.linalg.norm(inv(a)-f)<1e-9
    return state(plot('Linear approximation',series(t,f,'Input'),series(t,linear,f'M = {M}'),xlabel='t',ylabel='Amplitude'),plot('Nonlinear approximation',series(t,f,'Input'),series(t,nonlin,f'M = {M}'),xlabel='t',ylabel='Amplitude'),plot('Coefficient magnitudes',series(np.arange(len(a)),abs(a),'Original',stems=True),series(np.arange(len(a)),abs(s),'Retained',points=True),xlabel='Index',ylabel='Magnitude'),metrics=[metric('Linear SNR',snr(f,linear),'dB'),metric('Nonlinear SNR',snr(f,nonlin),'dB')],note='Both approximations use exactly M real orthonormal coefficients. Nonlinear approximation selects the largest magnitudes, including deterministic tie-breaking.')


def approximation2d(view,p):
    f=input_image(p['image']);M=p['M'];kind=p['basis'];g,b=approximation_image(p['image'],kind,M);linear,_=approximation_image(p['image'],kind,M,False)
    if view in ['approx-speed','approx-effi']:
        counts=np.array([16,32,64,128,256,512,1024,2048,4096,8192,12000]);curves=[]
        for basis in ['Fourier','DCT','haar','db4']:
            a,_,_=transform_image(p['image'],basis);e=np.cumsum(np.sort(abs(a).ravel()**2)[::-1]);err=np.maximum((np.sum(a*a)-e[counts-1])/np.sum(f*f),1e-16);curves.append(series(counts,err,basis))
        return state(plot('Relative squared approximation error',*curves,xlabel='M real coefficients',ylabel='Relative squared error',logx=True,logy=True),image(g,f'{kind}, M = {M}'),metrics=[metric('SNR',snr(f,g),'dB')])
    return state(image(f,'Original'),image(linear,'Linear approximation'),image(g,'Largest-coefficient approximation'),metrics=[metric('M / N',M/f.size),metric('Linear SNR',snr(f,linear),'dB'),metric('Nonlinear SNR',snr(f,g),'dB')],note='Equal budgets count real orthonormal coefficients in every basis. Images are 128 × 128; periodic boundaries are used for Fourier and wavelets.')


@functools.lru_cache(maxsize=30)
def mesh_case(angular,variation):
    from anisotropic_models import fitted_mesh,simple_cartoon,ellipse_radius,interpolate_mesh,mesh_checks
    pts,tri=fitted_mesh(angular);v=simple_cartoon(*pts.T);region=ellipse_radius(*pts.T)<1;base=np.where(region,.64,.32);v=base+variation*(v-base)
    g=interpolate_mesh(pts,tri,v,128);y,x=np.mgrid[:128,:128]/127;r=ellipse_radius(x,y)<1;f=np.where(r,.64,.32)+variation*(simple_cartoon(x,y)-np.where(r,.64,.32));checks=mesh_checks(pts,tri)
    return pts,tri,f,g,checks

def mesh(view,p):
    pts,tri,f,g,c=mesh_case(p['angular'],p['variation']);lines=[]
    for a,b,d in tri:lines.append(arr(pts[[a,b,d,a]]))
    meshpanel={'kind':'mesh','title':'Contour-aligned triangulation','lines':lines,'extent':[0,1,1,0],'background':image(g,'')['src']}
    metrics=[metric('Triangles',len(tri)),metric('Maximum tangent angle',c['jump_triangle_max_tangent_angle_degrees'],'°'),metric('SNR',snr(f,g),'dB')]
    if view=='triangulation-anisoiso':
        zoom={**meshpanel,'title':'Boundary detail','extent':[.67,.86,.66,.40]}
        return state(meshpanel,zoom,metrics=metrics,note='The detail shows the same mesh and continuous P1 interpolant near the right side of the elliptical jump. Long edges follow the tangent; the narrow direction resolves the jump.')
    if view=='triangulation-peppers':
        import struct
        from anisotropic_models import simple_cartoon,ellipse_radius,interpolate_mesh
        v=simple_cartoon(*pts.T);base=np.where(ellipse_radius(*pts.T)<1,.64,.32);v=base+p['variation']*(v-base)
        coords=np.round(pts*65535).astype('<u2');values=np.round(np.clip(v,0,1)*255).astype('uint8');cells=tri.astype('<u2')
        payload=b'P1M1'+struct.pack('<HHH',128,len(pts),len(tri))+coords.tobytes()+values.tobytes()+cells.tobytes()
        decoded=interpolate_mesh(coords/65535,tri,values/255,128);jpeg,stream,_=C.encode_image_at_budget(f,len(payload),'JPEG2000');budget=max(len(payload),len(stream))
        return state(image(f,'Original'),meshpanel,image(decoded,'Decoded P1 image'),image(jpeg,'JPEG-2000'),metrics=[metric('Equal transmitted budget',8*budget/f.size,'bpp'),metric('P1 SNR',snr(f,decoded),'dB'),metric('JPEG-2000 SNR',snr(f,jpeg),'dB')],note='The mesh stream explicitly contains uint16 positions, uint8 nodal intensities, uint16 triangle indices, and a header. Both streams are padded to the same transmitted byte budget. No optimal mesh entropy coding is assumed.')
    return state(image(f,'Piecewise-smooth input'),meshpanel,image(g,'Continuous P1 interpolant'),metrics=metrics,note='Thin triangles follow the singularity tangent. Smooth variations remain on both sides; the middle panel overlays the actual P1 interpolant, not the original image.')


def directional(view,p):
    from curvelet_frame import windows,analyze,synthesize,denoise
    if view=='curvelet-denoising':
        f=input_image(p['image'],64);sig=.08;y=f+sig*noise(f.shape);bank,_=windows(64);g=denoise(y,sig,p['factor'],bank);w=C.stationary_denoise(y,p['factor']*sig,'hard',level=2)
        assert np.linalg.norm(synthesize(analyze(f,bank),bank)-f)<1e-9
        return state(image(y,'Noisy input'),image(w,'Translation-invariant wavelets'),image(g,'Periodic curvelet frame'),metrics=[metric('Wavelet SNR',snr(f,w),'dB'),metric('Curvelet SNR',snr(f,g),'dB')],note='The documented oversampled Parseval frame is used, not the fast wrapping transform. Both thresholds use the same noise level; curvelet thresholds account for each atom’s norm.')
    if view=='curvelets-vs-wavelets':
        r=p['scale'];offset=np.deg2rad(p['angle']);theta=np.linspace(0,2*np.pi,501);panels=[]
        for directional_ in [False,True]:
            curves=[series(np.cos(theta),np.sin(theta),'Singularity boundary',color='#203746')]
            for k,a in enumerate(np.linspace(0,2*np.pi,16,endpoint=False)):
                t=np.linspace(0,2*np.pi,50);angle=a+np.pi/2+offset if directional_ else np.pi/2*round((a+np.pi/2)/(np.pi/2));long=1/np.sqrt(r) if directional_ else 1/r;short=1/r
                u=long*np.cos(t);v=short*np.sin(t);x=np.cos(a)+u*np.cos(angle)-v*np.sin(angle);y=np.sin(a)+u*np.sin(angle)+v*np.cos(angle)
                curves.append(series(x,y,'Localization' if k==0 else '',color='#146c78'))
            panels.append(plot('Parabolic curvelets' if directional_ else 'Isotropic wavelets',*curves,xlabel='x',ylabel='y',equal=True,xlim=[-1.6,1.6],ylim=[-1.6,1.6]))
        return state(*panels,note='Ellipses indicate localization, not exact compact supports. Curvelet length scales as 1/√r and width as 1/r. The orientation control offsets each elongated atom from the tangent, making a mismatch visible.')
    n=128;bank,info=windows(n);radius=p['scale'];angle=np.deg2rad(p['angle']);idx=min(range(1,len(info)),key=lambda i:abs(np.log2(info[i]['radius']/radius))*10+abs(np.angle(np.exp(1j*(info[i]['angle']-angle)))))
    win=bank[idx];atom=fft.fftshift(fft.ifft2(win,norm='ortho').real)
    if view=='curvelet-discretization':
        a=np.deg2rad(p['angle']);u=np.arange(-3,4)/np.sqrt(radius);v=np.arange(-6,7)/radius;U,V=np.meshgrid(u,v);X=U*np.cos(a)-V*np.sin(a);Y=U*np.sin(a)+V*np.cos(a)
        return state(plot('Parabolic translation grid',series(X.ravel(),Y.ravel(),'Positions',points=True),xlabel='x',ylabel='y',equal=True),image(atom,'Selected directional atom',signed=True),note='Tangential spacing scales as 1/√r; normal spacing as 1/r. The lattice rotates with the atom orientation.')
    return state(image(atom,'Spatial atom (real part)',signed=True),image(fft.fftshift(win),'Frequency window'),metrics=[metric('Radial scale',info[idx]['radius']),metric('Frequency angle',np.rad2deg(info[idx]['angle']),'°')],note='A smooth annular/angular window gives a parabolically scaled directional atom. Spatial oscillation is normal to its elongated envelope; displayed windows come from the same Parseval frame used for denoising.')


def compression(view,p):
    f=input_image(p['image']);budget=int(p['rate']*f.size/8)
    if view=='compression':
        T=p['rate'];a,s=C.wavelet_array(f,'db2',3);q=np.sign(a)*np.floor(abs(a)/T);b=np.sign(q)*(abs(q)+.5)*T;g=C.wavelet_inverse(b,s,'db2');M=np.count_nonzero(q)
        return state(image(f,'Original'),waveimage((q!=0).astype(float),s,'Quantized support'),image(g,'Decoded image'),metrics=[metric('Nonzero coefficients',M),metric('SNR',snr(f,g),'dB'),metric('Step T',T)],note='An explicit support codec: dead-zone quantization followed by midpoint reconstruction. The rate control changes the quantization step; the encoded support size is reported, without claiming JPEG-2000 rates.')
    j,bj,_=C.encode_image_at_budget(f,budget,'JPEG');w,bw,_=C.encode_image_at_budget(f,budget,'JPEG2000')
    return state(image(f,'Original'),image(j,'JPEG'),image(w,'JPEG-2000'),metrics=[metric('JPEG actual rate',len(bj)*8/f.size,'bpp'),metric('JPEG-2000 actual rate',len(bw)*8/f.size,'bpp'),metric('JPEG SNR',snr(f,j),'dB'),metric('JPEG-2000 SNR',snr(f,w),'dB')],note='Real encode/decode operations, including file headers. Each codec targets the same byte budget; its achieved rate is shown because not every target is attainable.')

FUNCTIONS={k:globals()[k] for k in ['acquisition','quantization','fourier','coding','wavelet1d','wavelet2d','models','approximation1d','approximation2d','mesh','directional','compression']}
