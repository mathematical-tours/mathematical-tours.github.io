from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


f=flower(256);M=1024;fig,a=canvas(1,3,height=3.4)
for ax,g,name in zip(a.ravel(),[f,approx(f,M,linear=True),approx(f,M)],['Original','Fixed coarse coefficients','Largest coefficients']):image(ax,g,name+(f'\nM = {M}; {snr(f,g):.1f} dB' if name!='Original' else ''))
finish(fig,__file__,'Linear and nonlinear approximations use the same orthonormal Daubechies-4 basis and exactly 1024 of 65536 real coefficients. The linear mask is the fixed 32 by 32 coarse block; nonlinear selection sorts magnitudes with a deterministic tie rule.',data_sources=['data/flower.png'],parameters={'size':256,'M':M,'wavelet':'db4','levels':4},checks={'nonlinear_error_no_larger':bool(np.linalg.norm(f-approx(f,M))<=np.linalg.norm(f-approx(f,M,linear=True)))})
