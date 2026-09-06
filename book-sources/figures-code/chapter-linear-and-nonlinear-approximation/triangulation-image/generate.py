from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from anisotropic_models import *
p,tri=fitted_mesh();v=simple_cartoon(*p.T);g=interpolate_mesh(p,tri,v);y,x=np.mgrid[:256,:256]/255;f=simple_cartoon(x,y)
fig,a=canvas(1,3,height=3.6,width=11.5)
image(a[0,0],f,'Smooth regions separated by a jump')
mesh_plot(a[0,1],p,tri,'Thin triangles along the boundary')
image(a[0,2],g,f'Continuous P1 interpolation\n{snr(f,g):.1f} dB')
checks=mesh_checks(p,tri);checks['P1_snr_db']=snr(f,g)
finish(fig,__file__,'The cartoon has different smooth nonaffine intensity functions inside and outside an elliptical jump. Nearly equilateral triangles fill the smooth regions; close contour-following layers create thin triangles whose long edges follow the tangent. Explicit strip connectivity prevents edge flips across the contour. The longest edges of all 128 jump-crossing triangles deviate by less than 5.4 degrees from the local tangent, with aspect ratios above 13.5. The third panel is the actual continuous P1 nodal interpolant. The contour is known analytically in this controlled example.',parameters={'angular_samples':64,'normal_radial_layers':[.91,.955,.985,.997,1.003,1.015,1.045,1.09]},checks=checks,arrays={'points':p,'triangles':tri,'values':v,'interpolation':g,'original':f})
