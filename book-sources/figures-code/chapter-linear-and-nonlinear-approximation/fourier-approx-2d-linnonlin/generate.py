from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


f=flower(128);fig,a=canvas(2,3,height=5.8)
for row,linear in enumerate([True,False]):
 for ax,M in zip(a[row],[256,1024,4096]):image(ax,approx(f,M,'Fourier',linear),f'{"Linear" if linear else "Nonlinear"}: M = {M}')
finish(fig,__file__,'The top row keeps a fixed disk of low spatial frequencies in a real Fourier basis. The bottom row chooses the largest M real coefficients of the same transform. Each column has an identical coefficient budget, image, and display range.',data_sources=['data/flower.png'],parameters={'size':128,'M':[256,1024,4096]})
