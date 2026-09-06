from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *

n=256;y,x=np.mgrid[:n,:n]/n
square=((abs(x-.5)<.25)&(abs(y-.5)<.25)).astype(float)
fig,axs=canvas(2,3,height=6.6);errors=[]
for col,(name,f) in enumerate([('Square',square),('Flower',flower(n)),('Flower detail',flower_detail(n))]):
    a,s=wavelet_array(f,'db2',4)
    image(axs[0,col],f,name);wavelet_image(axs[1,col],a,s,'Wavelet coefficients')
    errors.append(float(np.max(abs(wavelet_inverse(a,s,'db2')-f))))
assert max(errors)<1e-12
finish(fig,__file__,'Retains the original image-over-coefficients arrangement. A synthetic square, the supplied flower image and a central flower detail compare simple edges with photographic structure. Every lower panel is computed from the image directly above, with the same transform and display scale.',parameters={'wavelet':'db2','levels':4},checks={'inverse_max_errors':errors})
