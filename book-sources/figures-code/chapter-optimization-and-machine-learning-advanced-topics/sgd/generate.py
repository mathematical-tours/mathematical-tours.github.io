from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


from optimization_models import optimization_runs
epoch,runs,checks=optimization_runs();fig,a=canvas(height=4,width=7);ax=a[0,0];ax.plot(epoch,runs['SGD'],label='SGD');ax.plot(epoch,runs['Batch'],'--',color=INK,label='Batch gradient descent');ax.set(yscale='log',xlabel='Equivalent passes through the data',ylabel='Objective gap',title='Logistic objective: equal gradient-evaluation budgets');ax.legend();ax.grid(True)
finish(fig,__file__,'One full-gradient update costs n individual gradients and is compared with n sampled SGD updates per horizontal unit. Both optimize exactly the same quadratically regularized logistic objective, with gaps measured against one high-accuracy numerical reference. SGD uses independent uniform sampling and tau=(1/L_i)/(1+k/(3n)).',parameters={'seed':20260906,'epochs':60,'gap_floor':1e-16},checks=checks,arrays={'epochs':epoch,'SGD':np.array(runs['SGD']),'Batch':np.array(runs['Batch'])})
