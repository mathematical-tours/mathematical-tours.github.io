from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


f=flower();fig,a=canvas(1,3,height=3.6)
for ax,M in zip(a.ravel(),[512,2048,8192]):g=approx(f,M);image(ax,g,f'M = {M}\nSNR = {snr(f,g):.1f} dB')
finish(fig,__file__,'Largest-coefficient db4 reconstructions at increasing budgets use the same acquired flower photograph. SNR is measured on the unclipped reconstruction; display uses the common unit intensity range.',data_sources=['data/flower.png'],parameters={'size':256,'M':[512,2048,8192],'wavelet':'db4','levels':4})
