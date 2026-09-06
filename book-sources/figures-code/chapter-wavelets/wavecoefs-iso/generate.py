from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *

f=flower(256);w='db2';L=4
a=np.apply_along_axis(lambda v:np.concatenate(pywt.wavedec(v,w,level=L,mode='periodization')),1,f)
a=np.apply_along_axis(lambda v:np.concatenate(pywt.wavedec(v,w,level=L,mode='periodization')),0,a)
b,s=wavelet_array(f,w,L);fig,axs=canvas(1,2,height=4.3)
tensor_wavelet_image(axs[0,0],a,L,'Anisotropic tensor product')
wavelet_image(axs[0,1],b,s,'Isotropic wavelet pyramid',labels=True)
assert abs(np.sum(a*a)-np.sum(f*f))<1e-8
finish(fig,__file__,'Restores the original two packed coefficient layouts using the supplied flower image. Full horizontal and vertical separators identify independent tensor-product scales on the left; nested quadrant boundaries identify recursive low-low subdivision on the right. The filter family and display map are identical.',parameters={'wavelet':w,'levels':L},checks={'tensor_energy_error':float(abs(np.sum(a*a)-np.sum(f*f)))})
