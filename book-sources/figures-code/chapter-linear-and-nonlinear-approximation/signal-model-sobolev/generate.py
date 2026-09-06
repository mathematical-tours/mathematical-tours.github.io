from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
n=256;noise=np.random.default_rng(20260906).standard_normal((n,n));fig,a=canvas(1,4,height=3.1,width=11);energies=[];arrays={'noise':noise}
for ax,sigma in zip(a.ravel(),[1.5,3,6,12]):
    f=ndimage.gaussian_filter(noise,sigma,mode='wrap');f=.5+.15*(f-f.mean())/f.std();energies.append(.5*np.sum(gradient(f)**2));image(ax,f,rf'$\sigma={sigma:g}$ pixels');arrays[str(sigma)]=f
assert np.all(np.diff(energies)<0)
finish(fig,__file__,'One seeded pure white-noise image is progressively convolved with wider periodic Gaussian kernels. For comparison, each result is affinely normalized to the same mean and standard deviation. The measured gradient energy decreases with smoothing even after this contrast normalization.',parameters={'size':n,'seed':20260906,'sigma_pixels':[1.5,3,6,12],'display_mean':.5,'display_std':.15},checks={'gradient_energies':energies},arrays=arrays)
