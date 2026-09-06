from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
f=flower(256);d=gradient(f);dy,dx=d;scale=float(np.quantile(np.abs(d),.985));gamma=.45
def tone(z):return np.sign(z)*np.minimum(np.abs(z)/scale,1)**gamma
rgb=np.stack((.5+.5*tone(dx),np.full_like(dx,.5),.5+.5*tone(dy)),axis=-1)
fig,a=canvas(1,3,height=3.7,width=10);image(a[0,0],f,r'$f$');image(a[0,1],rgb,r'$\nabla f$: red $D_xf$, blue $D_yf$');obj=image(a[0,2],np.linalg.norm(d,axis=0),r'$|\nabla f|$',vmin=0,vmax=.3);fig.colorbar(obj,ax=a[0,2],shrink=.7)
rng=np.random.default_rng(5);p=rng.standard_normal(d.shape);err=abs(np.sum(d*p)+np.sum(f*divergence(p)));assert err<1e-10
finish(fig,__file__,'The middle RGB image encodes signed horizontal differences in red and signed vertical differences in blue. A shared contrast curve sign(z) min(|z|/s,1)^0.45, with s equal to the pooled 98.5th absolute-gradient percentile, strengthens weak components and clips only the display tails. Red and blue channels equal 0.5 plus half this signed curve; green stays at 0.5. The gradient arrays themselves are unchanged. Zero gradient is neutral gray; increasing/decreasing channels distinguish signs. The final panel is the actual Euclidean magnitude. Periodic div=-gradient-adjoint is verified.',parameters={'signed_channel_scale':float(scale),'scale_percentile':98.5,'display_gamma':gamma},checks={'negative_adjoint_error':float(err),'display_saturated_component_fraction':float(np.mean(abs(d)>scale))},arrays={'dx':dx,'dy':dy,'gradient_rgb':rgb})
