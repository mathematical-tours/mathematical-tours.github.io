from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=flower(256,True);fig,axs=canvas(1,3,height=3.6,width=10)
image(axs[0,0],f,'RGB image: 3 measured channels')
ax=axs[0,1];ax.set_axis_off();ax.set_title('Multispectral image: 32 channels (schematic)')
y,x=np.mgrid[:64,:64]/64
wavelength=np.linspace(400,1000,32)
cube=np.array([.2+.5*np.exp(-((l-470)/100)**2)*(x>.45)+.3*np.exp(-((l-780)/130)**2)*((x-.4)**2+(y-.5)**2<.2**2) for l in wavelength])
for j in range(31,-1,-1):
    ax.imshow(cube[j],cmap='gray',vmin=0,vmax=1,extent=(j*.014,1+j*.014,j*.012,1+j*.012),zorder=32-j)
    ax.plot([j*.014,1+j*.014],[1+j*.012]*2,color=MUTED,lw=.3,zorder=33-j)
ax.set(xlim=(-.02,1.6),ylim=(-.22,1.6));ax.text(.65,-.17,r'$f(x_1,x_2,\lambda_\ell),\quad \ell=1,\ldots,32$',ha='center')
axs[0,2].plot(wavelength,cube[:,32,40],'-o',ms=3);axs[0,2].set(title='Synthetic spectrum at one pixel',xlabel='Wavelength (nm)',ylabel='Intensity',ylim=(0,1));axs[0,2].grid(True)
finish(fig,__file__,'The RGB panel is the supplied flower photograph. The 32-band cube is a schematic synthetic spectral field with two wavelength-dependent components; it is not inferred from RGB channels and does not claim real multispectral measurements. This separates the spatial dimension d=2 from the channel count s.',data_sources=['data/flower.png'],parameters={'spectral_channels':32,'wavelength_nm':[400,1000],'synthetic_spectrum':True},arrays={'synthetic_spectral_cube':cube,'wavelength_nm':wavelength})
