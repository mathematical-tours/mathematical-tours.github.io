from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from applications_data import polygon,plot_polygon
z=polygon();rng=np.random.default_rng(20260908)
w=z+.025*(rng.standard_normal(len(z))+1j*rng.standard_normal(len(z)))
w=w.mean()+np.exp(.2j*np.pi)*(w-w.mean())
descriptor=lambda x:abs(np.fft.fft(x)[2:])/abs(np.fft.fft(x)[1])
assert np.allclose(descriptor(z),descriptor((2+3j)*z+4))
fig,axes=plt.subplots(1,2,figsize=(8.7,3.8))
plot_polygon(axes[0],z,color=TEAL,label='original');plot_polygon(axes[0],w,color=RED,ls='--',label='perturbed and rotated')
axes[0].set(aspect='equal',title='Polygons');axes[0].legend(frameon=False,fontsize=8);axes_box(axes[0],'real coordinate','imaginary coordinate')
for x,color,marker,label in [(z,TEAL,'o','original'),(w,RED,'x','modified')]:axes[1].plot(np.arange(2,len(z)),descriptor(x),marker,color=color,ms=4,label=label)
axes_box(axes[1],'frequency index $k$',r'$|\widehat\Pi[k]|/|\widehat\Pi[1]|$');axes[1].set_title('Magnitude descriptors');axes[1].legend(frameon=False,fontsize=8)
save(fig,'chapter-fourier-applications','fourier-descriptors')
