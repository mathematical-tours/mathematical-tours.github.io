"""Periodic, undecimated parabolic curvelet frame with an exact squared partition.

This transparent reference construction is not the fast wrapping implementation:
translations are deliberately oversampled on the full pixel grid. Smooth compact
annular/angular windows have radial width O(r) and transverse width O(sqrt(r)).
The low-frequency family and pointwise normalization make synthesis Parseval.
"""
from common import *

def bump(t):
    z=np.zeros_like(t);ok=abs(t)<1;z[ok]=np.exp(1-1/(1-t[ok]**2));return z

def windows(n=128):
    k=fft.fftfreq(n)*n;y,x=np.meshgrid(k,k,indexing='ij');r=np.hypot(x,y);theta=np.arctan2(y,x);bank=[bump(r/8)];info=[{'radius':0,'angle':0}]
    for radius in 2.**np.arange(2,int(np.log2(n))):
        angles=int(8*2**np.ceil(np.log2(radius/4)/2));radial=bump(np.log2(np.maximum(r,1e-20)/radius))
        for angle in np.arange(angles)*2*np.pi/angles:
            delta=np.angle(np.exp(1j*(theta-angle)));bank.append(radial*bump(delta/(2*np.pi/angles)));info.append({'radius':float(radius),'angle':float(angle)})
    total=np.sqrt(np.sum(np.array(bank)**2,axis=0));assert total.min()>0
    return [w/total for w in bank],info

def analyze(f,bank):return [fft.ifft2(fft.fft2(f,norm='ortho')*w,norm='ortho') for w in bank]
def synthesize(c,bank):return fft.ifft2(sum(fft.fft2(a,norm='ortho')*w for a,w in zip(c,bank)),norm='ortho').real

def denoise(y,sigma,factor,bank):
    c=analyze(y,bank)
    for a,w in zip(c[1:],bank[1:]):a*=abs(a)>factor*sigma*np.sqrt(np.mean(w*w))
    return synthesize(c,bank)
