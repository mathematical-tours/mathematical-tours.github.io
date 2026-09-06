from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from scipy.spatial.distance import cdist

rng=np.random.default_rng(20260906)
x=rng.uniform(-2,2,(75,2))
truth=lambda z:np.sin(1.8*z[:,0])*np.cos(1.3*z[:,1])+.25*z[:,0]
y=truth(x)+.12*rng.standard_normal(len(x))
t=np.linspace(-2,2,140);X,Y=np.meshgrid(t,t);grid=np.c_[X.ravel(),Y.ravel()]
widths=[.1,.5,1,5];ridge=.03;fig,axs=canvas(1,5,height=3.4,width=13)
vmin,vmax=-1.5,1.5
obj=axs[0,0].scatter(x[:,0],x[:,1],c=y,cmap='coolwarm',vmin=vmin,vmax=vmax,s=22,edgecolors=INK,linewidths=.25)
axs[0,0].set(title='Observed samples')
residuals=[]
for ax,width in zip(axs[0,1:],widths):
    K=np.exp(-cdist(x,x,'sqeuclidean')/(2*width**2))
    coef=np.linalg.solve(K+ridge*np.eye(len(x)),y)
    fitted=(np.exp(-cdist(grid,x,'sqeuclidean')/(2*width**2))@coef).reshape(X.shape)
    ax.imshow(fitted,origin='lower',extent=[-2,2,-2,2],cmap='coolwarm',vmin=vmin,vmax=vmax,interpolation='bilinear')
    ax.set_title(rf'$\sigma={width:g}$')
    residuals.append(float(np.linalg.norm((K+ridge*np.eye(len(x)))@coef-y)))
for ax in axs.ravel():ax.set(xlim=(-2,2),ylim=(-2,2),aspect='equal',xlabel=r'$x_1$',ylabel=r'$x_2$',xticks=[-2,0,2],yticks=[-2,0,2])
fig.colorbar(obj,ax=axs.ravel().tolist(),shrink=.7,label='Response / fitted value',pad=.015)
assert max(residuals)<1e-10
finish(fig,__file__,'Restores the original five-panel two-dimensional example: observed colored samples, then four Gaussian-kernel fits with the original bandwidths 0.1, 0.5, 1 and 5. All panels use the same coordinates and response color scale. Coefficients solve (K+0.03 I)c=y for one seeded noisy sample, and predictions evaluate the same kernel expansion on the grid.',parameters={'seed':20260906,'samples':75,'input_dimension':2,'kernel_widths':widths,'ridge_matrix_diagonal':ridge,'color_limits':[vmin,vmax]},checks={'linear_system_residuals':residuals})
