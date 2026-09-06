from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from scipy.optimize import minimize_scalar,brentq
fig,a=canvas(1,3,height=3.6,width=11);xx,yy=np.meshgrid(np.linspace(-1.3,2.2,240),np.linspace(-1.3,1.8,240));F=lambda x,y:(x-1.8)**2+2*(y-.3)**2
E=F(xx,yy);theta=minimize_scalar(lambda t:F(np.cos(t),np.sin(t)),bounds=(0,np.pi/2),method='bounded').x
p2=np.array([np.cos(theta),np.sin(theta)]);p1=np.array([1.,0.]);t=np.linspace(0,2*np.pi,801)
for ax,q,p,title in [(a[0,0],2,p2,'Smooth boundary'),(a[0,1],1,p1,'Corner contact')]:
 r=(abs(np.cos(t))**q+abs(np.sin(t))**q)**(-1/q);ax.fill(r*np.cos(t),r*np.sin(t),color=TEAL,alpha=.14);ax.plot(r*np.cos(t),r*np.sin(t),color=TEAL)
 level=F(*p);ax.contour(xx,yy,E,levels=sorted([level,level+.45,level+1.1,level+1.9]),colors=MUTED,linewidths=.7);ax.plot(*p,'o',color=RED,ms=6);ax.annotate(r'$x^\star$',p,xytext=p+[-.3,-.32],arrowprops={'arrowstyle':'-','color':RED},color=RED);ax.set(title=title,aspect='equal',xlim=(-1.15,2.1),ylim=(-1.1,1.5),xticks=[-1,0,1,2],yticks=[-1,0,1])
# n.x=1: the minimum l2 and l1 balls touch at different computed points.
ax=a[0,2];normal=np.array([1.,.6]);c=1.;u=normal/(normal@normal);v=np.array([1.,0.]);x=np.linspace(-1.2,2.1,300);ax.plot(x,(c-normal[0]*x)/normal[1],color=INK,label=r'$x_1+0.6x_2=1$')
r=np.linalg.norm(u);ax.plot(r*np.cos(t),r*np.sin(t),color=TEAL,label=r'Minimal $\ell^2$ ball');rho=1/(abs(np.cos(t))+abs(np.sin(t)));ax.plot(rho*np.cos(t),rho*np.sin(t),color=GOLD,label=r'Minimal $\ell^1$ ball')
for p,label in [(u,r'$x_2^\star$'),(v,r'$x_1^\star$')]:ax.plot(*p,'o',color=RED,ms=6);ax.annotate(label,p,xytext=p+[-.1,.16],color=RED)
ax.set(title='Affine constraint',aspect='equal',xlim=(-1.15,2.1),ylim=(-1.1,1.5),xticks=[-1,0,1,2],yticks=[-1,0,1]);ax.legend(loc='lower left',fontsize=8)
assert abs(normal@u-1)<1e-12 and abs(normal@v-1)<1e-12
finish(fig,__file__,'The first two panels show exact constraint contact points for F(x)=(x1-1.8)^2+2(x2-.3)^2. At the l1 corner, the supporting outward normal lies in the diamond normal cone. For x1+0.6x2=1, the minimal Euclidean ball touches at (1,.6)/1.36 and the minimal l1 ball at (1,0); both contacts and their correctly scaled balls are displayed.',checks={'l2_contact':p2.tolist(),'l1_contact':p1.tolist(),'affine_l2_contact':u.tolist(),'affine_l1_contact':v.tolist()})
