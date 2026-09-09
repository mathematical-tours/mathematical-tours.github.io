from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
n=48;x=np.linspace(-3,3,n,endpoint=False);X,Y=np.meshgrid(x,x)
f=3*(1-X)**2*np.exp(-X**2-(Y+1)**2)-10*(X/5-X**3-Y**5)*np.exp(-X**2-Y**2)-np.exp(-(X+1)**2-Y**2)/3
f+=.55*np.random.default_rng(20260908).standard_normal((n,n))
k=np.arange(-n//2,n//2);K,L=np.meshgrid(k,k);g=np.exp(-(K*K+L*L)/(2*1.8**2));g/=g.sum()
smoothed=np.fft.ifft2(np.fft.fft2(f)*np.fft.fft2(np.fft.ifftshift(g))).real
assert np.isclose(f.sum(),smoothed.sum())
fig=plt.figure(figsize=(11.5,3.4))
for panel,(Z,title) in enumerate([(f,'Noisy surface $f$'),(g,'Unit-sum Gaussian $g$'),(smoothed,r'Circular convolution $f*g$')]):
    ax=fig.add_subplot(1,3,panel+1,projection='3d')
    ax.plot_surface(X,Y,Z,cmap='viridis',rstride=1,cstride=1,linewidth=0,antialiased=True,rasterized=True)
    ax.set(xlabel='$x$',ylabel='$y$',title=title,xticks=[-3,0,3],yticks=[-3,0,3])
    if panel!=1:ax.set_zlim(-7,9)
    ax.view_init(28,-58);ax.set_box_aspect((1,1,.68));ax.tick_params(pad=0,labelsize=8)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False
save(fig,'chapter-discrete-fourier-transform','two-dimensional-convolution')
