from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from anisotropic_models import *
p,tri=fitted_mesh();fig,a=canvas(1,2,height=4.7,width=10)
mesh_plot(a[0,0],p,tri,'Uniform in smooth regions; thin along the edge',background=True)
x0,x1,y0,y1=.33,.57,.72,.96
a[0,0].add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fill=False,edgecolor=GOLD,lw=1.5))
mesh_plot(a[0,1],p,tri,'Magnified boundary region',background=True)
a[0,1].set(xlim=(x0,x1),ylim=(y1,y0))
for spine in a[0,1].spines.values():spine.set_edgecolor(GOLD);spine.set_linewidth(1.5)
# Show the tangent and normal at a point on the actual contour, in equal units.
theta=1.80;point=.5+ELLIPSE_AXES*np.array([np.cos(theta),np.sin(theta)])
tangent=ELLIPSE_AXES*np.array([-np.sin(theta),np.cos(theta)]);tangent/=np.linalg.norm(tangent)
normal=(point-.5)/ELLIPSE_AXES**2;normal/=np.linalg.norm(normal)
a[0,1].annotate('',xy=point+.045*tangent,xytext=point-.045*tangent,arrowprops={'arrowstyle':'<->','color':GOLD,'lw':2})
a[0,1].text(point[0]-.065,point[1]+.045,'long: tangent',color=INK,fontsize=11)
finish(fig,__file__,'The same mesh as Figure 5.20 is shown over its smoothly varying nodal interpolant. The orange box is magnified with equal horizontal and vertical units. The long triangle edges track the red singularity curve; their short altitudes cross it. The orange arrow indicates the tangent direction. Away from the contour, the mesh is nearly equilateral, without radial spokes.',parameters={'zoom_box':[x0,x1,y0,y1],'equal_axis_units':True},checks=mesh_checks(p,tri),arrays={'points':p,'triangles':tri})
