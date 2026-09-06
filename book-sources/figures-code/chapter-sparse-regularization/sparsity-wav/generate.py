from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


f=flower();c,s=wavelet_array(f);fig,a=canvas(1,3,height=3.6);image(a[0,0],f,'Photograph');wavelet_image(a[0,1],c,s,'Wavelet coefficients',gain=30);d=np.sort(abs(c).ravel())[::-1];curve(a[0,2],np.arange(1,len(d)+1),d,'Ordered magnitudes','Rank','Magnitude');a[0,2].set(xscale='log',yscale='log')
finish(fig,__file__,'The coefficient image and magnitude decay are both computed from the displayed flower. Logarithmic axes reveal many small coefficients without implying that they are exactly zero.',data_sources=['data/flower.png'],parameters={'size':256,'wavelet':'db4','levels':4})
