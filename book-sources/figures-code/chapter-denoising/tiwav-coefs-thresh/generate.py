from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case(n=512);cs=pywt.swtn(y,'db4',level=1,trim_approx=True,norm=False);a=cs[1]['da'];b=hard(a,2.7*.08);limit=4*.08;fig,axs=canvas(1,2,height=4)
image(axs[0,0],a,r'Level $j=-8$, orientation $H$',signed=True,cmap='gray',vmin=-limit,vmax=limit);image(axs[0,1],b,r'Hard threshold: $T=2.7\sigma$',signed=True,cmap='gray',vmin=-limit,vmax=limit)
finish(fig,__file__,'Uses a 512 by 512 flower image so the first detail level is j=-8 in the chapter’s J=-log2(N) convention. The H band is high-pass in the first array coordinate. The same signed color range [-4 sigma,4 sigma] is used, with saturation only for display, before and after exact hard thresholding; equality at the threshold maps to zero.',data_sources=['data/flower.png'],parameters={'size':512,'seed':20260906,'sigma':.08,'wavelet':'db4','j':-8,'threshold_sigma':2.7},checks={'retained_count':int(np.count_nonzero(b))})
