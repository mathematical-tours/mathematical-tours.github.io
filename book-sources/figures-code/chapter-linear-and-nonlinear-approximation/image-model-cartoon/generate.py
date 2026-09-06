from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
n=256;y,x=np.mgrid[:n,:n]/n;rng=np.random.default_rng(20260906);f=np.full((n,n),.5);fig,a=canvas(1,4,height=3.1,width=11);counts=[8,32,128,512];arrays={};energies=[]
for k in range(1,max(counts)+1):
    cx,cy=rng.uniform(-.1,1.1,2);r=rng.uniform(.04,.24);ratio=rng.uniform(.6,1.5);theta=rng.uniform(0,np.pi)
    u=(x-cx)*np.cos(theta)+(y-cy)*np.sin(theta);v=-(x-cx)*np.sin(theta)+(y-cy)*np.cos(theta)
    f[(u/r)**2+(v/(ratio*r))**2<=1]=rng.uniform(.08,.92)
    if k in counts:
        image(a[0,counts.index(k)],f.copy(),f'{k} opaque leaves');arrays[str(k)]=f.copy();energies.append(tv(f)/n)
finish(fig,__file__,'A progressive dead-leaves model overlays opaque ellipses with independently sampled positions, sizes, orientations and gray values. Later leaves occlude earlier ones. Every snapshot is piecewise constant with a finite union of smooth boundary arcs; intersections and occlusions create corners. Total variation is measured, not assumed monotone under occlusion.',parameters={'size':n,'seed':20260906,'leaf_counts':counts,'radius_range':[.04,.24]},checks={'tv_per_grid_side':energies},arrays=arrays)
