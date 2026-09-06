from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


f=cartoon();c,s=wavelet_array(f);M=1024;b=keep_largest(c,M);g=wavelet_inverse(b,s);fig,a=canvas(1,3,height=3.5);image(a[0,0],f,'Cartoon image');image(a[0,1],g,f'M = {M} retained coefficients');wavelet_image(a[0,2],c,s,'Wavelet coefficients',gain=30)
finish(fig,__file__,'A computed cartoon approximation and its actual orthogonal wavelet coefficients show concentration around edge curves. The signed reconstruction is not clipped before SNR measurement, and the coefficient display is explicitly logarithmic.',parameters={'size':256,'M':M,'wavelet':'db4','levels':4},checks={'snr_db':snr(f,g)})
