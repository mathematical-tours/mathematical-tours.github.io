"""Contour-aligned triangulations of a cartoon with smooth, nonaffine regions."""
from common import *
from scipy.spatial import Delaunay
import matplotlib.tri as mtri

ELLIPSE_CENTER=np.array([.5,.5])
ELLIPSE_AXES=np.array([.28,.34])


def ellipse_radius(x,y):
    return np.sqrt(((x-.5)/.28)**2+((y-.5)/.34)**2)


def simple_cartoon(x,y):
    """Two smooth functions, separated by a nonzero elliptical jump."""
    exterior=.20+.16*x+.10*y+.11*np.sin(2*np.pi*x)*np.cos(2*np.pi*y)
    interior=.59+.12*np.sin(2*np.pi*x+.4)+.10*np.cos(2*np.pi*y-.3)+.10*x*y
    return np.where(ellipse_radius(x,y)<1,interior,exterior)


def fitted_mesh(angular=64):
    # A quasiuniform triangular lattice fills the smooth regions. Only the
    # contour neighborhood is refined in the normal direction; no radial fan
    # connects the boundary to the image center or to the square boundary.
    theta=np.linspace(0,2*np.pi,angular,endpoint=False)
    unit=np.c_[np.cos(theta),np.sin(theta)]
    layers=np.array([.91,.955,.985,.997,1.003,1.015,1.045,1.09])
    rings=np.vstack([ELLIPSE_CENTER+r*unit*ELLIPSE_AXES for r in layers])
    lattice=[];spacing=1/13
    for row,y in enumerate(np.arange(spacing/2,1,spacing*np.sqrt(3)/2)):
        for x in np.arange(spacing/2+(row%2)*spacing/2,1,spacing):
            r=ellipse_radius(x,y)
            if r<.84 or r>1.16:lattice.append((x,y))
    t=np.linspace(0,1,17)
    boundary=np.vstack([np.c_[t,np.zeros_like(t)],np.c_[t,np.ones_like(t)],np.c_[np.zeros_like(t[1:-1]),t[1:-1]],np.c_[np.ones_like(t[1:-1]),t[1:-1]]])
    points=np.vstack([rings,lattice,boundary]);tri=[]
    lattice_indices=np.arange(len(rings),len(rings)+len(lattice))
    radii=ellipse_radius(*np.asarray(lattice).T)
    inner=np.r_[np.arange(angular),lattice_indices[radii<.84]]
    tri.extend(inner[Delaunay(points[inner]).simplices])
    # Triangulate the exterior independently, removing the convex polygonal
    # hole. Its ring edges are separated from the quasiuniform lattice.
    outer=np.r_[np.arange(7*angular,8*angular),lattice_indices[radii>1.16],np.arange(len(rings)+len(lattice),len(points))]
    exterior=Delaunay(points[outer]).simplices
    tri.extend(outer[exterior[np.any(exterior>=angular,axis=1)]])
    # Explicit connectivity is essential here: unconstrained Delaunay flips
    # would replace some contour edges in the very thin transition strip.
    for j in range(len(layers)-1):
        for k in range(angular):
            h=(k+1)%angular;a=j*angular;b=a+angular
            tri.extend([[a+k,b+k,b+h],[a+k,b+h,a+h]])
    triangles=np.asarray(tri,dtype=int)
    e=points[triangles[:,1]]-points[triangles[:,0]]
    f=points[triangles[:,2]]-points[triangles[:,0]]
    flip=e[:,0]*f[:,1]-e[:,1]*f[:,0]<0
    triangles[flip]=triangles[flip][:,[0,2,1]]
    # The two polygonal chains bracketing the jump must remain mesh edges.
    edges={tuple(sorted(e)) for tri in triangles for e in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]}
    for layer in (3,4):
        for k in range(angular):assert tuple(sorted((layer*angular+k,layer*angular+(k+1)%angular))) in edges
    checks=mesh_checks(points,triangles)
    assert checks['jump_triangle_max_tangent_angle_degrees']<8
    assert checks['jump_triangle_min_aspect_ratio']>12
    assert checks['smooth_triangle_median_aspect_ratio']<2.2
    return points,triangles


def mesh_checks(points,triangles):
    vertices=points[triangles];centers=vertices.mean(axis=1)
    e=vertices[:,1]-vertices[:,0];f=vertices[:,2]-vertices[:,0]
    twice_area=e[:,0]*f[:,1]-e[:,1]*f[:,0]
    assert np.all(twice_area>0)
    assert abs(twice_area.sum()/2-1)<1e-12
    edges=np.stack([vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,1],vertices[:,0]-vertices[:,2]],axis=1)
    length2=np.sum(edges**2,axis=-1);longest=edges[np.arange(len(edges)),np.argmax(length2,axis=1)]
    # Longest edge divided by its altitude: a dimensionless shape measure.
    aspect=length2.max(axis=1)/twice_area
    radii=ellipse_radius(vertices[:,:,0],vertices[:,:,1]);crossing=(radii.min(axis=1)<1)&(radii.max(axis=1)>1)
    normal=(centers-ELLIPSE_CENTER)/(ELLIPSE_AXES**2)
    tangent=np.c_[-normal[:,1],normal[:,0]]
    cosine=np.abs(np.sum(longest*tangent,axis=1))/(np.linalg.norm(longest,axis=1)*np.linalg.norm(tangent,axis=1))
    angles=np.degrees(np.arccos(np.clip(cosine,0,1)))
    smooth=(radii.max(axis=1)<.84)|(radii.min(axis=1)>1.16)
    return {'vertices':len(points),'triangles':len(triangles),'domain_area':float(twice_area.sum()/2),
            'jump_triangles':int(crossing.sum()),'jump_triangle_max_tangent_angle_degrees':float(angles[crossing].max()),
            'jump_triangle_median_aspect_ratio':float(np.median(aspect[crossing])),
            'jump_triangle_min_aspect_ratio':float(aspect[crossing].min()),
            'smooth_triangle_median_aspect_ratio':float(np.median(aspect[smooth]))}


def interpolate_mesh(points,triangles,values,n=256):
    y,x=np.mgrid[:n,:n]/(n-1)
    interpolator=mtri.LinearTriInterpolator(mtri.Triangulation(*points.T,triangles),values)
    result=interpolator(x,y)
    assert not np.any(np.ma.getmaskarray(result))
    return np.asarray(result)


def mesh_plot(ax,p,tri,title='',background=False):
    if background:ax.tripcolor(*p.T,tri,simple_cartoon(*p.T),shading='gouraud',cmap='gray',vmin=0,vmax=1,alpha=.28,rasterized=True)
    ax.triplot(*p.T,tri,color=TEAL,lw=.43)
    theta=np.linspace(0,2*np.pi,1000)
    ax.plot(.5+.28*np.cos(theta),.5+.34*np.sin(theta),color=RED,lw=1.15)
    ax.set(aspect='equal',xlim=(0,1),ylim=(1,0),title=title,xticks=[],yticks=[])
