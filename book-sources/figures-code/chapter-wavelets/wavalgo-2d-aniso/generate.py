from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *

f=flower(256);L=4;w='db2'
a=np.apply_along_axis(lambda v:np.concatenate(pywt.wavedec(v,w,level=L,mode='periodization')),1,f)
b=np.apply_along_axis(lambda v:np.concatenate(pywt.wavedec(v,w,level=L,mode='periodization')),0,a)
fig,axs=canvas(1,3,height=3.7)
image(axs[0,0],f,'Image $f$')
tensor_wavelet_image(axs[0,1],a,L,'Transform each row',axes=(1,))
tensor_wavelet_image(axs[0,2],b,L,'Transform each column')
assert abs(np.sum(b*b)-np.sum(f*f))<1e-8
finish(fig,__file__,'Restores the flower photograph and the original row-then-column presentation. Separators follow the four-level one-dimensional packing after the first pass and its full tensor-product grid after the second pass. Every coefficient is computed, not copied from the old coefficient image.',parameters={'wavelet':w,'levels_per_axis':L},checks={'energy_preserved':True})
