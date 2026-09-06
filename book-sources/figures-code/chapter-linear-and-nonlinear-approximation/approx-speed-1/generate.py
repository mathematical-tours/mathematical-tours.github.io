from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


fig,a=canvas(1,4,height=3.4)
for ax,(name,f) in zip(a.ravel(),models().items()):image(ax,f,name)
finish(fig,__file__,'These four precisely defined test images are reused unchanged in the next approximation-error figure: smooth, cartoon, the supplied flower and the author-requested mandrill.',data_sources=['data/flower.png'],parameters={'size':256})
