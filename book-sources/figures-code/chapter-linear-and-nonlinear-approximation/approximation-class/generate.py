from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import models
rng=np.random.default_rng(20260906);noise=rng.standard_normal((128,128));bv=tv_denoise(noise,.6,steps=700);bv=.5+.16*bv/max(bv.std(),1e-12)
d=models();panels={'Smooth':d['Smooth'],'Bounded variation (ROF)':bv,'Cartoon':d['Cartoon'],'Natural image':d['Flower']};fig,a=canvas(1,4,height=3.1,width=11)
for ax,(name,f) in zip(a.ravel(),panels.items()):image(ax,f,name)
finish(fig,__file__,'The BV panel is obtained by solving the ROF total-variation denoising problem on a seeded pure white-noise image. An affine rescaling makes its piecewise-flat structure visible. The classes overlap: smooth and cartoon functions also belong to BV.',parameters={'noise_side':128,'seed':20260906,'ROF_lambda':.6,'ROF_iterations':700,'display_mean':.5,'display_std':.16},checks={'noise_tv':tv(noise),'ROF_tv':tv(bv)},arrays={'noise':noise,'ROF_image':bv})
