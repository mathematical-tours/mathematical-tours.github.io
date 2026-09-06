from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


fig,a=canvas(1,2,height=4)
for ax,reverse in zip(a.ravel(),[False,True]):
 diagram_axis(ax)
 if not reverse:
  for y,j in [(.75,1),(.3,2)]:box(ax,.16,y,rf'$\dot x_{{\ell_{j}}}$',width=.22);arrow(ax,(.29,y),(.69,.53),rf'$J_{{k,\ell_{j}}}$')
  box(ax,.81,.53,r'$\dot x_k$',width=.2);ax.text(.5,.07,r'$\dot x_k=\sum_\ell J_{k\ell}\dot x_\ell$',ha='center');ax.set_title('Forward: sum contributions of parents')
 else:
  box(ax,.16,.53,r'$\bar x_k$',width=.2)
  for y,j in [(.75,1),(.3,2)]:box(ax,.81,y,rf'$\bar x_{{m_{j}}}$',width=.22);arrow(ax,(.68,y),(.28,.53),rf'$J_{{m_{j},k}}^{{\mathsf{{T}}}}$')
  ax.text(.5,.07,r'$\bar x_k=\sum_m J_{mk}^{\mathsf{T}}\bar x_m$',ha='center');ax.set_title('Reverse: sum contributions of children')
finish(fig,__file__,'Forward arrows transport tangent vectors with local Jacobians; reverse arrows transport adjoint vectors with transposed local Jacobians. Sums include every parent or child, respectively. The order of matrix multiplication matches the vector-valued chain rule.')
