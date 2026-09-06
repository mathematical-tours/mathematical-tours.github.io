from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
t=np.linspace(0,1,2400);phi,psi,z=pywt.Wavelet('db4').wavefun(level=10)
fig,axs=canvas(3,2,height=6.2)
shifts_by_level=[]
for row,k in enumerate([2,4,8]):
    curve(axs[row,0],t,np.sqrt(2)*np.cos(2*np.pi*k*t),rf'Fourier atom: $k={k}$',r'$t$')
    j=row+4;shifts=[0,4,8] if j==4 else ([1,10,19] if j==5 else [2,18,34,50]);shifts_by_level.append(shifts)
    for n in shifts:
        v=2**(j/2)*np.interp(2**j*t-n,z,psi,left=0,right=0);axs[row,1].plot(t,v,label=rf'$n={n}$',lw=1.2)
    axs[row,1].set(title=rf'Translated Daubechies-4 atoms, $j={j}$',xlabel=r'$t$',xlim=(0,1));axs[row,1].legend(ncol=2,fontsize=8);axs[row,1].grid(True)
finish(fig,__file__,'Each right panel overlays several integer translations of the same actual Daubechies-4 wavelet at a fixed scale. The displayed atoms are 2^(j/2) psi(2^j t-n); supports stay inside the common unit interval. Distinct colors identify translations, while successive rows change frequency or scale.',parameters={'wavelet':'db4','frequencies':[2,4,8],'wavelet_levels':[4,5,6],'translations':shifts_by_level})
