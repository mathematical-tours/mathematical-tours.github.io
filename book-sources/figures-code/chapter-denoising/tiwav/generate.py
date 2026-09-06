from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=flower(256);cs=pywt.swtn(f,'db2',level=2,trim_approx=True,norm=False);fig,axs=canvas(2,4,height=5.8);image(axs[0,0],f,'Image');image(axs[1,0],cs[0]/4,'Coarsest approximation / 4')
for row,detail in enumerate(cs[1:]):
    for col,(key,label) in enumerate([('da','H'),('ad','V'),('dd','D')],1):
        image(axs[row,col],detail[key],f'Level {2-row}, {label}',signed=True,cmap='gray')
finish(fig,__file__,'Shows the actual undecimated two-dimensional Daubechies-2 transform. Every detail array has the original image shape; no subsampling is performed. H means high-pass in the first array coordinate (key da), V in the second (ad), and D in both (dd), matching the tensor-factor convention. Detail grays are signed; the coarse array is divided by four for display.',data_sources=['data/flower.png'],parameters={'size':256,'wavelet':'db2','levels':2,'stationary_normalization':False})
