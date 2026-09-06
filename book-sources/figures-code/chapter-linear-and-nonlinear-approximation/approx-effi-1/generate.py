from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import transform
f=flower(256);M=4096;fig,a=canvas(1,3,height=3.5,width=10);records={}
for ax,kind in zip(a.ravel(),['Fourier','DCT','Wavelet']):
    c,inv=transform(f,kind);g=inv(keep_largest(c,M));e=np.linalg.norm(f-g)**2/np.linalg.norm(f)**2;image(ax,g,f'{kind}: M = {M}\nRelative squared error {e:.4f}');records[kind]=float(e)
    assert abs(np.sum((f-g)**2)-np.sum((c-keep_largest(c,M))**2))<1e-8
finish(fig,__file__,'All three orthonormal representations retain exactly 4096 real coefficients from the same 256 by 256 flower image. Fourier counts individual real sine/cosine atoms; wavelets use four periodic db4 levels. At this larger common budget, wavelets have the lowest measured error. This observation concerns this image and budget, not a universal ordering of bases.',parameters={'size':256,'M':M,'wavelet':'db4','levels':4},checks={'relative_squared_errors':records,'discarded_energy_identity':True})
