from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=flower(256);fig,axs=canvas(1,4,height=3.2)
levels=[2,3,4,16]
for ax,K in zip(axs.ravel(),levels):
    q=np.clip(np.floor((K-1)*f+.5),0,K-1)/(K-1)
    image(ax,q,f'{K} gray levels')
    assert len(np.unique(q))<=K
finish(fig,__file__,'Nearest-level quantization uses K equally spaced reproduction values between zero and one: Q(u)=clip(floor((K-1)u+1/2),0,K-1)/(K-1). All panels use the same flower image and display range. The finite endpoints and tie convention are explicit.',data_sources=['data/flower.png'],parameters={'levels':levels},checks={'number_of_levels_verified':True})
