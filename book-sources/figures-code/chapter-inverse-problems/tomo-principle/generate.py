from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *

f=phantom(400);theta=np.pi/6;normal=np.array([np.cos(theta),np.sin(theta)]);tangent=np.array([-normal[1],normal[0]]);offset=.18
s=np.linspace(-2,2,600);line=offset*normal[:,None]+tangent[:,None]*s
positions=np.linspace(-1,1,1000);projection=phantom_projection(positions,theta);value=float(phantom_projection(np.array([offset]),theta)[0])
fig,a=canvas(1,2,height=4,width=10);ax=a[0,0];ax.imshow(f,cmap='gray',vmin=0,vmax=1,extent=[-1,1,-1,1]);ax.plot(*line,color=RED,lw=1.4);ax.annotate('',xy=.8*normal,xytext=(0,0),arrowprops={'arrowstyle':'->','color':GOLD});ax.text(*(.86*normal),r'$n_\theta$',color=GOLD);ax.set(xlim=(-1,1),ylim=(-1,1),aspect='equal',title='Modified Shepp–Logan phantom',xlabel=r'$x_1$',ylabel=r'$x_2$')
ax=a[0,1];curve(ax,positions,projection,r'One projection: $\theta=\pi/6$',r'$t$',r'$\mathcal{R}f(t,\theta)$');ax.vlines(offset,0,value,color=RED,ls='--');ax.plot(offset,value,'o',color=RED);ax.set_xlim(-1,1)
assert abs(normal@tangent)<1e-12 and np.max(abs(normal@line-offset))<1e-12
# Integrating a Radon profile gives the total mass of the ten ellipses.
mass=np.sum(np.pi*PHANTOM_ELLIPSES[:,0]*PHANTOM_ELLIPSES[:,1]*PHANTOM_ELLIPSES[:,2]);error=abs(np.trapezoid(projection,positions)-mass);assert error<.002
finish(fig,__file__,'The conventional modified Shepp-Logan test object is generated from ten ellipses. The red line has normal n_theta and signed offset t; the red dot is its actual line integral on the exact analytic one-dimensional Radon projection. Integrating the projection recovers the phantom mass. Ellipse conventions follow the standard modified Shepp-Logan model described in https://www.mathworks.com/help/images/ref/phantom.html .',parameters={'theta':float(theta),'signed_offset':offset,'ellipses':PHANTOM_ELLIPSES.tolist()},checks={'projection_mass_error':float(error),'marked_projection_value':value})
