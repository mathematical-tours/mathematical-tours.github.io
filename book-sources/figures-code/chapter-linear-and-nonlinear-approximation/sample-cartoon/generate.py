from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
y,x=np.mgrid[:256,:256]/256;g=.2+.6*(y>.5+.18*np.sin(2*np.pi*x));f=cartoon();fig,a=canvas(1,3,height=3.2,width=11)
image(a[0,0],cartoon_frame(),'Feline Follies (1919)')
image(a[0,1],g,'Smooth graph discontinuity');image(a[0,2],ndimage.gaussian_filter(f,2,mode='wrap'),'Smoothed mathematical cartoon')
finish(fig,__file__,'The first panel is a real grayscale animation still from Feline Follies (1919), attributed by Wikimedia Commons to Pat Sullivan and listed there as public domain. Source: https://commons.wikimedia.org/wiki/File:Felix_1919.jpg . The other panels illustrate a mathematically defined jump along a smooth graph and an actual Gaussian-smoothed cartoon. An animation drawing and the mathematical cartoon class are distinguished.',parameters={'blur_sigma_pixels':2})
