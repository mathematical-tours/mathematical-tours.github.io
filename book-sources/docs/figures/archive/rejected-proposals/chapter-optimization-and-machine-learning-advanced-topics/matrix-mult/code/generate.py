from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


fig,a=canvas(1,2,height=3.6)
for ax,reverse in zip(a.ravel(),[False,True]):
 diagram_axis(ax);ax.set_title('Forward accumulation' if not reverse else 'Reverse accumulation')
 ax.text(.5,.78,r'$A_3(A_2A_1)$' if not reverse else r'$(A_3A_2)A_1$',ha='center',fontsize=20)
 ax.text(.5,.54,r'$A_1:40\times100,\quad A_2:20\times40$'+'\n'+r'$A_3:1\times20$',ha='center',fontsize=12)
 ax.text(.5,.28,'82,000 scalar multiplications' if not reverse else '4,800 scalar multiplications',ha='center',color=TEAL,fontsize=13)
 ax.text(.5,.07,r'$n_0\sum_{k=1}^{t-1}n_kn_{k+1}$' if not reverse else r'$n_t\sum_{k=0}^{t-2}n_kn_{k+1}$',ha='center',fontsize=15)
finish(fig,__file__,'Both parenthesizations compute the same chain Jacobian A3 A2 A1. Exact dense multiplication counts for dimensions (100,40,20,1) illustrate the general formulas. Counts refer to scalar multiplications and exclude additions and the shared cost of evaluating local Jacobians.',parameters={'dimensions':[100,40,20,1]},checks={'forward_multiplications':82000,'reverse_multiplications':4800})
