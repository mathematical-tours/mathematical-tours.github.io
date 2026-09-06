from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=flower(256);fig,axs=canvas(1,4,height=3.4)
for ax,bins in zip(axs.ravel(),[128,64,32,16]):
    k=256//bins;a=f.reshape(bins,k,bins,k).mean(axis=(1,3));image(ax,np.repeat(np.repeat(a,k,axis=0),k,axis=1),rf'$j={-int(np.log2(bins))}$: {bins} × {bins}')
finish(fig,__file__,'Every displayed image is a true Haar orthogonal projection obtained by averaging the supplied flower image on dyadic squares. Increasing j enlarges the squares and reduces the number of degrees of freedom; all four panels share the same intensity range.',data_sources=['data/flower.png'],parameters={'square_grids':[128,64,32,16]})
