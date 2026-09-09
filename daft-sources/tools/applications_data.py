"""Shared inputs for related applications figures, with deterministic data."""
from pathlib import Path
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]


def flower(n=128):
    with Image.open(ROOT/'data/flower.png') as im:
        im=im.convert('L').resize((n,n),Image.Resampling.LANCZOS)
        return np.asarray(im,dtype=float)/255


def irregular_signal(n=128):
    rng=np.random.default_rng(20260908);m=n//4;x=np.linspace(0,1,m)
    return np.r_[1+.4*rng.random(m),np.sin(np.pi*x),1+rng.random(m),x*x]


def gaussian_1d(sigma, radius=None):
    if sigma==0:return np.array([0]),np.ones(1)
    radius=int(np.ceil(4*sigma)) if radius is None else radius
    k=np.arange(-radius,radius+1);g=np.exp(-k*k/(2*sigma*sigma))
    return k,g/g.sum()


def polygon():
    # Vertex coordinates from the original DAFT generator; orient counterclockwise.
    x=np.array([.6606,.5345,.3048,.2306,.1995,.2150,.3187,.4482,.7228,.7884,
                .7798,.6641,.4914,.4378,.4655,.5881,.7280,.7003,.5915])
    y=np.array([.0892,.0817,.1621,.2952,.5214,.7274,.8078,.8656,.8555,.7224,
                .5992,.6168,.6269,.5013,.3379,.2676,.1420,.0867,.0716])
    z=x+1j*y
    if np.sum((z.conj()*np.roll(z,-1)).imag)<0:z=z[::-1]
    return z


def polygon_evolution(z,multiplier,steps):
    F=np.fft.fft(z)
    return np.stack([np.fft.ifft(F*multiplier**m) for m in steps])


def plot_polygon(ax,z,**kwargs):
    closed=np.r_[z,z[0]]
    return ax.plot(closed.real,closed.imag,**kwargs)
