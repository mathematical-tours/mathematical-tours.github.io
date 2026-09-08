"""Geometry, random matrices, real-data learning, and optimization experiments."""
from numerics import *
from scipy.optimize import minimize,linprog
from sklearn.linear_model import LogisticRegression,lasso_path
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
from learning_models import iris,digits,pca as pca_fit,logistic_case


def controls(family,view):
    if family=='lq':return [slider('q','Exponent q',[.3,.5,.75,1,1.5,2,4,8]),slider('radius','Affine right-hand side' if view=='optim-geom' else 'Radius',[.5,.75,1,1.25,1.5])]
    if family=='certificates':return [slider('s','Support size s',[2,3,5,8,12]),slider('P','Measurements P',[16,24,32,40,48])]
    if family=='random':
        if view in ['marcenko-pastur','marcenko-pastur-laws']:return [slider('beta','Aspect ratio β = s / P',[.1,.2,.4,.6,.8]),slider('P','Rows P',[40,80,160,320])]
        if view in ['max-gaussian','max-gaussian-mean']:return [slider('N','Number of variables N',[8,16,32,64,128,256,512,1024])]
        return [slider('P','Measurements P',[12,20,28,36,44]),slider('s','Support size s',[1,2,4,6,8,10,12])]
    if family=='pca':return [slider('component','Second principal component',[2,3,4]),select('dataset','Dataset',['Iris','Digits'],1 if view in ['digits','digits-classes'] else 0)]
    if family=='clustering':return [slider('iteration','Lloyd iteration',[0,1,2,3,5,10,20,40]),slider('K','Clusters K',[2,3,4,6,8])]
    if family=='learning':
        if view=='hist-classif':return [slider('k','Neighbors R',[1,3,5,9,15,31,61]),select('dataset','Geometry',['Moons','Circles'])]
        if view=='losses':return [slider('margin','Negative margin s',[-3,-2,-1,0,.5,1,2,3])]
        if view=='separation-influ':return [slider('separation','Class separation',[.2,.5,1,1.5,2,3]),slider('C','Inverse regularization C',[.1,1,10,100])]
        return [slider('sigma','Kernel width σ',[.05,.1,.2,.4,.8,1.6]),slider('lam','Regularization λ',[.0001,.001,.01,.1,1])]
    if family=='stochastic':
        if view=='sgd-traject':return [slider('tau','Step size τ',[.02,.05,.1,.2,.5,1,1.5]),slider('sigma','Gradient noise σ',[.1,.3,.5,1])]
        return [slider('factor','Step-size multiplier',[.1,.3,.6,1,1.5,2]),slider('epochs','Data passes',[5,10,20,40,60])]
    raise KeyError(family)


def lq(view,p):
    q=p['q'];r=p['radius'];t=np.linspace(0,2*np.pi,801);x=np.cos(t);y=np.sin(t);rad=r*(abs(x)**q+abs(y)**q)**(-1/q)
    curves=[series(rad*x,rad*y,f'ℓq radius {r}')]
    if view=='optim-geom':
        u=np.linspace(-1.7,1.7,301);v=(r-u)/2;curves=[series(u,v,'x₁ + 2x₂ = b',dash=True)]
        for power,name in [(2,'ℓ²'),(1,'ℓ¹'),(q,'Selected q')]:
            a=np.array([1.,2.]);contact=r*a**(1/(power-1))/sum(a**(power/(power-1))) if power>1 else np.array([0,r/2])
            radius=np.sum(abs(contact)**power)**(1/power);rr=radius*(abs(x)**power+abs(y)**power)**(-1/power)
            curves.extend([series(rr*x,rr*y,name+' touching ball'),series([contact[0]],[contact[1]],name+' contact',points=True)])
        note='Each ball is the smallest scaled ball meeting the affine constraint. Contact points are computed analytically. For q ≤ 1 the minimizer lies on the coordinate with the largest constraint coefficient; q < 1 gives a nonconvex gauge.'
    else:note='The set is {x : |x₁|^q + |x₂|^q ≤ r^q}. It is convex exactly when q ≥ 1.'
    return state(plot('ℓq geometry',*curves,xlabel='x₁',ylabel='x₂',equal=True,xlim=[-1.8,1.8],ylim=[-1.8,1.8]),metrics=[metric('q',q),metric('Convex','yes' if q>=1 else 'no')],note=note)


@functools.lru_cache(maxsize=40)
def cert_case(s,P):
    rng=np.random.default_rng(SEED);full=rng.standard_normal((48,64));A=full[:P]/np.sqrt(P);I=np.arange(s)*3%64;sgn=np.where(np.arange(s)%2,1.,-1.);B=A[:,I];eta=A.T@B@np.linalg.solve(B.T@B,sgn);J=np.setdiff1d(np.arange(64),I);x=np.zeros(64);x[I]=sgn;y=A@x
    assert np.max(abs(eta[I]-sgn))<1e-8
    return A,x,y,I,J,eta

def certificates(view,p):
    s,P=p['s'],p['P'];A,x,y,I,J,eta=cert_case(s,P)
    if view=='homotopy':
        lammax=max(abs(A.T@y));lams=np.geomspace(lammax,.002*lammax,60);_,coefs,_=lasso_path(A,y,alphas=lams/P,max_iter=40000,tol=1e-9)
        return state(plot('Lasso solution path',*[series(lams,coefs[k],f'x{k+1}') for k in range(64) if np.max(abs(coefs[k]))>.08],xlabel='λ',ylabel='Coefficient',logx=True),metrics=[metric('Support size',s),metric('Measurements',P)],note='Each point solves ½‖Ax − y‖² + λ‖x‖₁. The solver’s mean-loss convention is compensated by dividing λ by P.')
    if view=='certif-convol':
        # Smooth continuous representation using a bandlimited convolution kernel.
        freq=np.arange(-P//2,P//2+1);loc=np.linspace(.12,.88,s);B=np.exp(-2j*np.pi*freq[:,None]*loc);weights=np.exp(-.5*(freq/(P/8))**2);B=weights[:,None]*B;signs=np.where(np.arange(s)%2,1.,-1.);v=B@np.linalg.solve(B.conj().T@B,signs);t=np.linspace(0,1,501);eta=(np.exp(2j*np.pi*t[:,None]*freq)*weights)@v
        return state(plot('Convolution precertificate',series(t,eta.real,'ηF'),series([0,1],[1,1],'+1',dash=True),series([0,1],[-1,-1],'−1',dash=True),series(loc,signs,'Support values',points=True),xlabel='Position t',ylabel='Certificate'),note='The Fuchs precertificate interpolates the support signs. The continuous extension uses the same weighted Fourier convolution operator; narrow spacing can violate the off-support bound.')
    bound=max(abs(eta[J]));return state(plot('Fuchs precertificate',series(np.arange(64),eta,'ηF',stems=True),series(I,np.sign(x[I]),'Support',points=True),series([0,63],[1,1],'+1',dash=True),series([0,63],[-1,-1],'−1',dash=True),xlabel='Index',ylabel='Certificate'),metrics=[metric('Off-support maximum',bound),metric('Nondegenerate',bound<1)],note='ηF = Aᵀ AI (AIᵀ AI)⁻¹ sign(xI). The same Gaussian matrix and nested support are used across the controls.')


@functools.lru_cache(maxsize=30)
def mp_case(beta,P):
    s=int(beta*P);B=np.random.default_rng(SEED).standard_normal((P,s))/np.sqrt(P);e=np.linalg.eigvalsh(B.T@B);return e,s/P

def random_case(view,p):
    if view in ['marcenko-pastur','marcenko-pastur-laws']:
        e,beta=mp_case(p['beta'],p['P']);lo=(1-np.sqrt(beta))**2;hi=(1+np.sqrt(beta))**2;x=np.linspace(lo,hi,401);d=np.sqrt(np.maximum((hi-x)*(x-lo),0))/(2*np.pi*beta*x);hist,edge=np.histogram(e,bins=20,range=(lo,hi),density=False);density=hist/(len(e)*np.diff(edge))
        return state(plot('Wishart eigenvalue distribution',series((edge[1:]+edge[:-1])/2,density,'Empirical',bars=True),series(x,d,'Marchenko–Pastur'),xlabel='Eigenvalue λ',ylabel='Density'),metrics=[metric('Actual s / P',beta),metric('Lower support edge',lo),metric('Upper support edge',hi),metric('Smallest eigenvalue',e[0])],note='B has independent N(0, 1/P) entries; the eigenvalues are those of BᵀB. All offered aspect ratios are below one, so the limiting law has no atom at zero.')
    if view in ['max-gaussian','max-gaussian-mean']:
        N=p['N'];z=np.random.default_rng(SEED).standard_normal((800,1024));vals=np.max(abs(z[:,:N]) if view=='max-gaussian-mean' else z[:,:N],axis=1);hist,edges=np.histogram(vals,35,density=True);grid=np.array([8,16,32,64,128,256,512,1024]);means=[np.max(abs(z[:,:n]) if view=='max-gaussian-mean' else z[:,:n],axis=1).mean() for n in grid]
        return state(plot('Distribution of the maximum',series((edges[:-1]+edges[1:])/2,hist,'800 trials',bars=True),xlabel='Maximum',ylabel='Density'),plot('Mean growth',series(grid,means,'Empirical mean'),series(grid,np.sqrt(2*np.log(grid)),'√(2 log N)',dash=True),xlabel='N',ylabel='Mean maximum',logx=True),metrics=[metric('Mean',vals.mean()),metric('Standard deviation',vals.std())],note=('Maximum absolute value' if view=='max-gaussian-mean' else 'Maximum (without absolute value)')+' of independent standard Gaussians; 800 trials with common random draws.')
    P,s=p['P'],p['s'];N=48;rng=np.random.default_rng(SEED);A=rng.standard_normal((P,N))/np.sqrt(P)
    if view=='rip-const':
        ks=np.arange(1,13);bounds=[]
        for k in ks:
            best=0
            for _ in range(80):
                I=rng.choice(N,k,replace=False);e=np.linalg.eigvalsh(A[:,I].T@A[:,I]);best=max(best,max(abs(e-1)))
            bounds.append(best)
        # Keep the maximum seen on every smaller support: these supports can be
        # embedded in larger ones, so the resulting lower bound is monotone.
        bounds=np.maximum.accumulate(bounds)
        return state(plot('Sampled RIP lower bounds',series(ks,bounds,'80 supports per size'),series([s],[bounds[s-1]],'Selected size',points=True),xlabel='Support size k',ylabel='Lower bound on δk'),note='Only sampled supports are tested. These are lower bounds on the worst-case RIP constants, never certified upper bounds.')
    counts=np.zeros(4);trials=32
    for _ in range(trials):
        A=rng.standard_normal((P,N));A/=np.linalg.norm(A,axis=0);I=rng.choice(N,s,replace=False);J=np.setdiff1d(np.arange(N),I);sgn=rng.choice([-1.,1.],s);B=A[:,I];G=B.T@B;cross=B.T@A[:,J];U=np.linalg.lstsq(G,cross,rcond=None)[0];den=1-np.linalg.norm(G-np.eye(s),1);weak=den>0 and np.linalg.norm(cross,1)<den;erc=np.linalg.norm(U,1)<1;fuchs=max(abs(sgn@U))<1;x=np.zeros(N);x[I]=sgn;sol=linprog(np.ones(2*N),A_eq=np.c_[A,-A],b_eq=A@x,bounds=(0,None),method='highs');ok=sol.success and np.linalg.norm(sol.x[:N]-sol.x[N:]-x)<1e-5;counts+=np.array([weak,erc,fuchs,ok])
    return state(plot('Empirical recovery criteria',series(np.arange(4),counts/trials,'Success fraction',bars=True),xlabel='0: weak ERC · 1: ERC · 2: Fuchs · 3: recovery',ylabel='Fraction',ylim=[0,1]),table('Monte Carlo counts',['Criterion','Successes / trials'],[[n,f'{int(c)} / {trials}'] for n,c in zip(['Weak ERC','ERC','Fuchs','Exact ℓ¹ recovery'],counts)]),note='32 independent Gaussian trials. Exact recovery is checked with linear programming. Finite trial counts explain irregular changes between neighboring parameters.')


@functools.lru_cache(maxsize=4)
def real_dataset(name):
    X,y=iris() if name=='Iris' else digits();m,s,V,Z=pca_fit(X);return X,y,m,s,V,Z

@functools.lru_cache(maxsize=2)
def trained_classifier(name):
    X,y=iris() if name=='Iris' else digits();train,test=train_test_split(np.arange(len(X)),test_size=.25,random_state=SEED,stratify=y);m,sv,V,Z=pca_fit(X[train]);model=LogisticRegression(C=3,max_iter=2500).fit(X[train],y[train]);return model,m,V,X[test],y[test],Z


def pca(view,p):
    X,y,m,s,V,Z=real_dataset(p['dataset']);comp=p['component']-1;curves=[]
    for label in np.unique(y):
        keep=np.flatnonzero(y==label)[::2 if len(y)>500 else 1];curves.append(series(Z[keep,0],Z[keep,comp],str(label),points=True))
    proj=plot('Principal-component projection',*curves,xlabel='PC1',ylabel='PC'+str(comp+1),equal=True)
    variance=s*s/(len(X)-1)
    if view=='cov':return state(image(np.cov(X,rowvar=False),'Empirical covariance',signed=True),plot('Covariance eigenvalues',series(np.arange(1,len(variance)+1),variance,'Eigenvalues',bars=True),xlabel='Principal component',ylabel='Variance'),proj)
    if view=='digits-classes':
        model,mean,axes,Xtest,ytest,trainZ=trained_classifier(p['dataset']);lo=trainZ[:,[0,comp]].min(axis=0)-.3;hi=trainZ[:,[0,comp]].max(axis=0)+.3
        gy,gx=np.mgrid[:72,:72]/71;grid=lo+(hi-lo)*np.c_[gx.ravel(),gy.ravel()];lift=mean+grid@axes[[0,comp]];prob=model.predict_proba(lift);assert np.max(abs(prob.sum(axis=1)-1))<1e-10
        palette=np.array([C.plt.get_cmap('tab10')(int(k))[:3] for k in model.classes_]);panel=image((prob@palette).reshape(72,72,3)[::-1],'Class-probability mixture on the PCA plane');panel['extent']=[float(lo[0]),float(hi[0]),float(hi[1]),float(lo[1])];panel['xlabel']='PC1';panel['ylabel']='PC'+str(comp+1)
        panel['colorLabels']=[{'label':str(k),'color':'#'+''.join(f'{round(255*v):02x}' for v in color)} for k,color in zip(model.classes_,palette)]
        return state(panel,table('Full-space classifier',['Quantity','Value'],[['Held-out accuracy',f'{model.score(Xtest,ytest):.3f}'],['Test points',len(ytest)]]),note='PCA and logistic regression are fitted on training data only. Each grid point is lifted from the selected two-PC plane into the full feature space before prediction. The color mixture uses all class probabilities; a projected test point need not have the same prediction as its affine-plane lift.')
    if view=='digits' and p['dataset']=='Digits':
        exemplars=[X[np.flatnonzero(y==k)[0]].reshape(8,8) for k in range(10)];sheet=np.vstack([np.hstack(exemplars[:5]),np.hstack(exemplars[5:])])
        return state(image(sheet,'Recorded handwritten digits 0–9'),proj,metrics=[metric('Images',len(X)),metric('Pixel features',X.shape[1])],note='These are real 8 × 8 handwritten digit images. The slider selects which principal component is paired with PC1; no labels are used to fit PCA.')
    return state(proj,plot('Explained variance',series(np.arange(1,len(s)+1),np.cumsum(s*s)/sum(s*s),'Cumulative fraction'),xlabel='Number of components',ylabel='Explained fraction',ylim=[0,1]),metrics=[metric('Observations',len(X)),metric('Dimensions',X.shape[1])],note='Real Iris measurements or handwritten digits from the repository. Data are centered; PCA is obtained by SVD, without per-feature variance normalization.')


@functools.lru_cache(maxsize=20)
def clustering_case(K,continuous):
    rng=np.random.default_rng(SEED)
    if continuous:
        pts=rng.random((1400,2));weights=.15+np.exp(-((pts[:,0]-.7)**2+(pts[:,1]-.35)**2)/.045);init=rng.random((K,2))
    else:
        X,y,*_=real_dataset('Iris');_,_,_,_,_,Z=real_dataset('Iris');pts=Z[:,:2];weights=np.ones(len(pts));init=pts[rng.choice(len(pts),K,replace=False)].copy()
    centers=init;history=[centers.copy()];energies=[]
    for _ in range(40):
        dist=np.sum((pts[:,None]-centers[None])**2,axis=2);lab=np.argmin(dist,axis=1);energies.append(float(np.sum(weights*dist[np.arange(len(pts)),lab])/weights.sum()));new=centers.copy()
        for k in range(K):
            I=lab==k
            if I.any():new[k]=np.average(pts[I],weights=weights[I],axis=0)
        centers=new;history.append(centers.copy())
    assert np.max(np.diff(energies))<1e-10
    return pts,weights,history,energies

def clustering(view,p):
    pts,weights,h,E=clustering_case(p['K'],view=='lloyd');k=p['iteration'];c=h[k];lab=np.argmin(np.sum((pts[:,None]-c[None])**2,axis=2),axis=1);curves=[series(pts[lab==i,0],pts[lab==i,1],f'Cluster {i+1}',points=True) for i in range(len(c))];curves.append(series(c[:,0],c[:,1],'Centroids',points=True,size=7,color='#182f3c'))
    return state(plot('Lloyd assignments and centroids',*curves,xlabel='x₁' if view=='lloyd' else 'PC1',ylabel='x₂' if view=='lloyd' else 'PC2',equal=True),plot('Quantization energy',series(np.arange(len(E)),E,'Energy'),series([min(k,39)],[E[min(k,39)]],'Current iteration',points=True),xlabel='Iteration',ylabel='Weighted squared distance'),metrics=[metric('Iteration',k),metric('Clusters',p['K'])],note='Weighted Lloyd iterations approximate a nonuniform continuous density using 1400 fixed quadrature samples.' if view=='lloyd' else 'Lloyd’s algorithm acts on the first two PCA coordinates of the real Iris measurements. Initialization stays fixed for each K.')


def classification_data(kind='Moons'):
    from sklearn.datasets import make_moons,make_circles
    return make_moons(220,noise=.17,random_state=SEED) if kind=='Moons' else make_circles(220,noise=.10,factor=.45,random_state=SEED)

def probability_panel(model,X,title,prob=True):
    lo=X.min(axis=0)-.25;hi=X.max(axis=0)+.25;y,x=np.mgrid[0:72,0:72]/71;grid=lo+(hi-lo)*np.c_[x.ravel(),y.ravel()];pr=model.predict_proba(grid)[:,1].reshape(72,72) if prob else model(grid).reshape(72,72)
    colors=np.array([[.18,.42,.62],[.75,.27,.24]]);rgb=(1-pr[...,None])*colors[0]+pr[...,None]*colors[1]
    panel=image(rgb[::-1],title);panel['probability']=True;panel['extent']=[float(lo[0]),float(hi[0]),float(hi[1]),float(lo[1])];panel['xlabel']='x₁';panel['ylabel']='x₂';return panel


def learning(view,p):
    if view=='losses':
        x=np.linspace(-3,3,401);s=p['margin'];loss=[('Binary',(x>=0).astype(float)),('Hinge',np.maximum(0,1+x)),('Logistic',np.logaddexp(0,x)),('Exponential',np.exp(x))]
        return state(plot('Loss as a function of negative margin',*[series(x,y,n) for n,y in loss],series([s,s],[0,5],'Current margin',dash=True),xlabel='s = −y ⟨x, β⟩',ylabel='Loss',ylim=[0,5]),table('Current losses',['Loss','Value'],[['Binary',int(s>=0)],['Hinge',f'{max(0,1+s):.3f}'],['Logistic',f'{np.logaddexp(0,s):.3f}'],['Exponential',f'{np.exp(s):.3f}']]),note='The unscaled logistic loss is log(1 + exp(s)); at zero it equals log 2. The binary loss counts zero margin as an error, as in this chapter.')
    if view=='hist-classif':
        X,y=classification_data(p['dataset']);model=KNeighborsClassifier(p['k']).fit(X,y);return state(plot('Training observations',*[series(X[y==i,0],X[y==i,1],f'Class {i}',points=True) for i in [0,1]],xlabel='x₁',ylabel='x₂',equal=True),probability_panel(model,X,'Neighbor vote probability P(y = 1 | x)'),metrics=[metric('Neighbors R',p['k'])],note='Class probability is the fraction of the R nearest training points in class 1. Every grid point is classified from the displayed fixed sample.')
    if view=='separation-influ':
        rng=np.random.default_rng(SEED);y=np.r_[np.zeros(100),np.ones(100)];X=rng.standard_normal((200,2))*.5;X[:,0]+=(2*y-1)*p['separation']/2;model=LogisticRegression(C=p['C']).fit(X,y)
        return state(plot('Training observations',*[series(X[y==i,0],X[y==i,1],f'Class {i}',points=True) for i in [0,1]],xlabel='x₁',ylabel='x₂',equal=True),probability_panel(model,X,'Logistic probability'),metrics=[metric('Training accuracy',model.score(X,y))],note='Separation moves the class centers while preserving every within-class displacement. C is the inverse strength of the quadratic penalty.')
    sigma=p['sigma'];lam=p['lam'];rng=np.random.default_rng(SEED)
    if view=='kernel':
        x=np.linspace(0,1,30);f=lambda t:np.sin(2*np.pi*t)+.35*np.cos(6*np.pi*t);y=f(x)+.15*rng.standard_normal(len(x));K=np.exp(-((x[:,None]-x[None,:])/sigma)**2/2);a=np.linalg.solve(K+lam*np.eye(len(x)),y);t=np.linspace(0,1,401);g=np.exp(-((t[:,None]-x)/sigma)**2/2)@a
        return state(plot('Gaussian-kernel ridge regression',series(t,f(t),'Ground truth',dash=True),series(x,y,'Observations',points=True),series(t,g,'Kernel estimate'),xlabel='x',ylabel='f(x)'),metrics=[metric('Kernel width',sigma),metric('λ',lam)],note='The coefficients solve (K + λI)a = y for the unnormalized squared-loss objective. Samples and observation noise remain fixed.')
    X,y=classification_data();K=np.exp(-np.sum((X[:,None]-X[None])**2,axis=2)/(2*sigma*sigma));model=LogisticRegression(C=1/lam,max_iter=2000).fit(K,y)
    class KernelModel:
        def predict_proba(self,g):return model.predict_proba(np.exp(-np.sum((g[:,None]-X[None])**2,axis=2)/(2*sigma*sigma)))
    return state(plot('Nonlinear classes',*[series(X[y==i,0],X[y==i,1],str(i),points=True) for i in [0,1]],xlabel='x₁',ylabel='x₂',equal=True),probability_panel(KernelModel(),X,'Gaussian-feature logistic probability'),note='Gaussian similarities to the training points are the features. Logistic fitting penalizes the feature coefficient vector; this finite feature model illustrates kernel nonlinear separation.')


@functools.lru_cache(maxsize=10)
def stochastic_case(factor):
    X,y,lam,loss,grad,grad_i=logistic_case();n,d=X.shape;fit=minimize(loss,np.zeros(d),jac=grad,method='L-BFGS-B',options={'gtol':1e-10,'ftol':1e-14});best=loss(fit.x);L=np.linalg.norm(X,2)**2/(4*n)+lam;Li=max(np.sum(X*X,axis=1))/4+lam;rng=np.random.default_rng(SEED);idx=rng.integers(n,size=60*n)
    states={k:np.zeros(d) for k in ['SGD','Average','SAG','Batch']};mean=np.zeros(d);table_=np.zeros((n,d));gg=np.zeros(d);runs={k:[loss(v)-best] for k,v in states.items()}
    for k,i in enumerate(idx):
        states['SGD']-=factor/Li/(1+k/(3*n))*grad_i(states['SGD'],i);mean+=(states['SGD']-mean)/(k+1)
        h=grad_i(states['SAG'],i);gg+=(h-table_[i])/n;table_[i]=h;states['SAG']-=factor*gg/(16*Li)
        if (k+1)%n==0:
            states['Batch']-=factor*grad(states['Batch'])/L
            for name,v in states.items():runs[name].append(max(loss(mean if name=='Average' else v)-best,1e-15))
    assert np.linalg.norm(grad(fit.x))<1e-5
    return runs

def stochastic(view,p):
    if view=='sgd-traject':
        tau=p['tau'];sig=p['sigma'];rng=np.random.default_rng(SEED);z=np.ones((120,200))*2
        for k in range(1,200):z[:,k]=z[:,k-1]-tau*(z[:,k-1]+sig*rng.standard_normal(120))
        variance=tau*sig*sig/(2-tau)
        return state(plot('Independent SGD trajectories',*[series(np.arange(200),z[i],f'Run {i+1}') for i in range(5)],xlabel='Iteration k',ylabel='xk'),plot('Empirical mean and standard deviation',series(np.arange(200),z.mean(axis=0),'Mean'),series(np.arange(200),z.std(axis=0),'Standard deviation'),series([0,199],[np.sqrt(variance)]*2,'Stationary standard deviation',dash=True),xlabel='Iteration k',ylabel='Value'),metrics=[metric('Stationary variance',variance)],note='For f(x)=x²/2, xk+1=(1−τ)xk−τσξk with independent standard Gaussian ξk. The stationary variance is τσ²/(2−τ), valid for 0<τ<2.')
    runs=stochastic_case(p['factor']);E=p['epochs'];keys=['Batch'] if view=='bgd' else (['SGD','Batch'] if view=='sgd' else ['SGD','Average','SAG','Batch'])
    return state(plot('Logistic objective gaps',*[series(np.arange(E+1),runs[k][:E+1],k) for k in keys],xlabel='Equivalent full data passes',ylabel='f(x) − f⋆',logy=True),metrics=[metric('Step multiplier',p['factor']),metric('Data passes',E)],note='All methods optimize the same regularized finite-sum logistic problem. A batch step and n individual gradients each cost one data pass. SGD has a decreasing step; averaging shares its iterates; SAG uses a stored gradient table.')

FUNCTIONS={'random':random_case,**{k:globals()[k] for k in ['lq','certificates','pca','clustering','learning','stochastic']}}
