from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


f=flower(128);mask=np.random.default_rng(20260906).random(f.shape)>.7;y=f*mask;x=y.copy();targets=[0,5,50,500];frames=[];energies=[]
for k in range(501):
 if k in targets:frames.append(x.copy());energies.append(float(.5*np.sum(gradient(x)**2)))
 x+=.24*laplacian(x);x[mask]=y[mask]
fig,a=canvas(1,4,height=3.5)
for ax,g,k in zip(a.ravel(),frames,targets):image(ax,g,f'Iteration {k}\n{snr(f,g):.1f} dB')
assert np.all(np.diff(energies)<=1e-8)
finish(fig,__file__,'Every update is the actual projected Sobolev step x <- projection_C(x+.24 Delta x). Projection restores all known pixels exactly, and the step is below 1/4. Reported energies decrease and measured-data consistency is checked.',data_sources=['data/flower.png'],parameters={'seed':20260906,'step':.24,'iterations':targets},checks={'sobolev_energies':energies,'known_pixel_error':float(max(np.max(abs(g[mask]-y[mask])) for g in frames))})
