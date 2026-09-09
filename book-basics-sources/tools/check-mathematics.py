#!/usr/bin/env python3
"""Independent checks of numerical examples and important corrected distinctions.
These checks support the editorial audit; they are not proofs of the theorems.
"""
from pathlib import Path
import itertools
import json
import math
import numpy as np
from PIL import Image
from scipy.optimize import linprog, linear_sum_assignment
from scipy.special import logsumexp

ROOT=Path(__file__).resolve().parents[1]
checks=[]
def check(name,condition):
    if not bool(condition): raise AssertionError(name)
    checks.append(name)

# Coding: check the actual finite example, Kraft mass, and bit strings.
message=[0,1,3,2,0,3,2,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,1,2,1]
codes=['001','01','1','000']
fixed=''.join(format(v,'02b') for v in message)
variable=''.join(codes[v] for v in message)
check('50-bit fixed code',fixed=='00011110001110100110100101101010010110101001011001')
check('42-bit prefix code',variable=='001010001001000110111010111101011110101101')
check('Empirical counts',np.array_equal(np.bincount(message),[2,9,12,2]))
check('Prefix condition',all(not b.startswith(a) for a,b in itertools.permutations(codes,2)))
p=np.bincount(message)/25
entropy=-sum(p*np.log2(p));length=sum(p*np.array(list(map(len,codes))))
check('Entropy and mean length',abs(entropy-1.6219011889093375)<1e-12 and abs(length-1.68)<1e-12)
check('Kraft sum',sum(2**(-len(c)) for c in codes)==1)
# A uniquely decodable non-prefix code: enumerate strings of varied lengths.
seen=set()
for n in range(11):
 for word in itertools.product(('0','01'),repeat=n):
  encoded=''.join(word);check_unique=encoded not in seen;seen.add(encoded)
  assert check_unique
checks.append('No collisions for {0,01} through ten source symbols (illustration only)')
for i,lens in enumerate(([2,3,4],[1,2,3])):
 ell=np.array(lens);p0=np.array([.2,.3,.5]);z=sum(2.**(-ell));q=2.**(-ell)/z
 kl=sum(p0*np.log2(p0/q));gap=sum(p0*ell)+sum(p0*np.log2(p0))
 check(f'Shannon KL identity, code {i+1}',abs(gap-(kl-math.log2(z)))<1e-12 and gap>=kl>=0)
check('Payload sizes',256*256*8==524288 and 256*256==65536 and 240*240/1024==56.25)
check('Binary expansion',int('10110011',2)==179)
patch=np.array([190,192,79,54,47,153,203,189,166])
check('Mean and median',patch.sum()==1273 and np.median(patch)==166)
check('Centered gradient example',(192-153)**2+(54-189)**2==19746)
check('Quantization bins',np.array_equal(np.arange(256)//64,np.repeat(np.arange(4),64)))
check('Additive primaries',np.array_equal(np.array([255,0,0])+[0,255,0],[255,255,0]))
a=np.linspace(0,255,500)
check('Power-map directions',np.all(a*a/255<=a+1e-12) and np.all(np.sqrt(255*a)>=a-1e-12))
check('Power-map inverses before rounding',np.allclose(np.sqrt(255*(a*a/255)),a))
folder=ROOT/'figures/images'
orig=np.asarray(Image.open(folder/'section10-original.png'))
check('Published transpose image',np.array_equal(orig.swapaxes(0,1),np.asarray(Image.open(folder/'section10-transpose.png'))))
check('Published clockwise rotation image',np.array_equal(np.rot90(orig,-1),np.asarray(Image.open(folder/'section10-rotation.png'))))
# Distinguish removal of regularization from sparse recovery.
A=np.array([[1.,1.]])
check('Large Lasso radius does not ensure recovery',np.allclose(A@np.array([1.,0.]),A@np.array([.5,.5])))
# A small, fixed compressed-sensing example, solved as a linear program.
rng=np.random.default_rng(38);N=64;P=36;M=4
A=rng.normal(size=(P,N))/np.sqrt(P);x0=np.zeros(N);x0[[2,18,31,49]]=[1.,-1.5,.7,2.]
r=linprog(np.ones(2*N),A_eq=np.c_[A,-A],b_eq=A@x0,bounds=(0,None),method='highs')
check('Basis pursuit recovers the sparse numerical example',r.success and np.max(abs(r.x[:N]-r.x[N:]-x0))<1e-8)
# Paris cost matrix, independently transcribed from the published illustration.
C=np.array([[12,10,31,27,10,30],[22,7,25,15,11,14],[19,7,19,10,15,15],
            [10,6,21,19,14,24],[15,23,14,24,31,34],[35,26,16,9,34,15]])
check('Paris example costs 65',C[np.arange(6),[4,1,5,0,2,3]].sum()==65)
check('Paris proposed optimum costs 64',C[np.arange(6),[4,5,1,0,2,3]].sum()==64)
check('Paris exhaustive optimum',min(sum(C[i,s[i]] for i in range(6)) for s in itertools.permutations(range(6)))==64)
u=np.array([2,3,3,0,5,2]);v=np.array([10,4,9,7,8,11])
check('Paris dual feasibility',np.all(u[:,None]+v[None,:]<=C))
check('Paris primal-dual equality',u.sum()+v.sum()==64)
check('Insertion-sort comparison count',70*69//2==2415)
for _ in range(100):
 x=np.sort(rng.normal(size=8));y=np.sort(rng.normal(size=8));cost=abs(x[:,None]-y[None,:])
 i,j=linear_sum_assignment(cost);assert np.isclose(cost[i,j].sum(),abs(x-y).sum())
checks.append('Sorted matching agrees with optimal assignment in 100 one-dimensional cases')
check('Squared W2 fails the triangle inequality',4>1+1)
# Explicit two-layer softmax-cross-entropy gradient versus finite differences.
d,m,k=3,4,2
x=np.array([.3,-.7,.5]);W=rng.normal(size=(m,d));b=np.array([2.,-2.,1.,-1.]);V=rng.normal(size=(k,m));c=rng.normal(size=k);y=np.array([0.,1.])
r=W@x+b;check('Gradient test avoids ReLU kinks',min(abs(r))>.1)
u=np.maximum(r,0);a=V@u+c;phat=np.exp(a-logsumexp(a));delta=phat-y;q=(V.T@delta)*(r>0)
analytic=np.r_[np.outer(q,x).ravel(),q,np.outer(delta,u).ravel(),delta]
weights=np.r_[W.ravel(),b,V.ravel(),c]
def loss(w):
 W=w[:m*d].reshape(m,d);b=w[m*d:m*d+m];V=w[m*d+m:m*d+m+k*m].reshape(k,m);c=w[-k:]
 a=V@np.maximum(W@x+b,0)+c
 return logsumexp(a)-a[1]
h=1e-6;fd=np.zeros_like(weights)
for j in range(len(weights)):
 wp=weights.copy();wm=weights.copy();wp[j]+=h;wm[j]-=h;fd[j]=(loss(wp)-loss(wm))/(2*h)
check('Backpropagation gradient versus centered differences',max(abs(fd-analytic))<2e-8)
check('Softmax translation invariance',np.allclose(np.exp(a-logsumexp(a)),np.exp(a+100-logsumexp(a+100))))
# The GAN identity is exact for finite strictly positive distributions.
p=np.array([.1,.3,.6]);q=np.array([.4,.4,.2]);mean=(p+q)/2;D=p/(p+q)
gan=sum(p*np.log(D)+q*np.log(1-D));js=.5*sum(p*np.log(p/mean))+.5*sum(q*np.log(q/mean))
check('GAN optimized discriminator / Jensen-Shannon identity',abs(gan-(-np.log(4)+2*js))<1e-12)
report={'checks_passed':len(checks),'checks':checks,'entropy_example':float(entropy),
        'backpropagation_max_error':float(max(abs(fd-analytic)))}
(ROOT/'reports/mathematics.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'{len(checks)} mathematical checks passed.')
