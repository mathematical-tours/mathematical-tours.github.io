from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


fig,a=canvas(height=4,width=9);ax=a[0,0];diagram_axis(ax)
for x,label in [(.07,r'$x$'),(.3,r'$W_1x$'),(.53,r'$\rho(W_1x)$'),(.76,r'$W_2^{\mathsf{T}}\rho(W_1x)$')]:box(ax,x,.37,label,width=.18,size=11)
box(ax,.95,.37,r'$+$',width=.06)
for left,right in [(.17,.2),(.4,.43),(.63,.66),(.86,.91)]:arrow(ax,(left,.37),(right,.37))
ax.plot([.07,.07,.95],[.46,.8,.8],color=RED,lw=1.2);arrow(ax,(.95,.8),(.95,.46),color=RED);ax.text(.48,.84,'Identity skip connection',ha='center',color=RED);ax.text(.5,.07,r'$h(x,\theta)=x+W_2^{\mathsf{T}}\rho(W_1x),\quad W_1,W_2\in\mathbb{R}^{q\times n}$',ha='center',fontsize=13)
finish(fig,__file__,'The residual block preserves the input dimension: W1 maps n to q, the activation acts in q dimensions, and W2 transpose maps q back to n before addition to x. The skip connection carries the original input without activation. Repeating the same block shares both weight matrices.')
