from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


from optimization_models import optimization_runs
epoch,runs,checks=optimization_runs();fig,a=canvas(height=4,width=7)
for name in ['SGD','SGA','SAG']:curve(a[0,0],epoch,np.log10(runs[name]),'Measured finite-sum optimization performance','Equivalent passes through the data',r'$\log_{10}(f(x_k)-f_{ref})$',label=name)
a[0,0].legend();finish(fig,__file__,'SGD, Cesàro-averaged SGD (SGA), and SAG use the same logistic problem and sampled-index stream. SGA averages its own auxiliary iterates. SAG maintains the exact mean of a table of past gradients, starts with zero entries, and updates the sampled entry before stepping; its direction is not called unbiased. Schedules and a conservative SAG step 1/(16 L_i) are explicit in the shared source.',parameters={'seed':20260906,'epochs':60,'gap_floor':1e-16},checks=checks,arrays={'epochs':epoch,**{k:np.array(runs[k]) for k in ['SGD','SGA','SAG']}})
