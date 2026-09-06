from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(height=3.2,width=11);ax=axs[0,0];diagram_axis(ax)
labels=['Input image','Wavelet\ntransform','Quantization','Code-block\nbit planes','Arithmetic\ncoding','Rate allocation\nand packets'];xs=[.08,.25,.42,.59,.76,.93]
for x,label in zip(xs,labels):box(ax,x,.64,label,width=.13,height=.19,size=10)
for a,b in zip(xs,xs[1:]):arrow(ax,(a+.075,.64),(b-.075,.64))
ax.text(.25,.36,'9/7 irreversible\n5/3 reversible',ha='center',fontsize=11,color=TEAL)
ax.text(.59,.36,'significance\nsign · refinement',ha='center',fontsize=11,color=TEAL)
ax.text(.93,.36,'quality layers\nprogressive decoding',ha='center',fontsize=11,color=TEAL)
ax.text(.5,.04,'Symmetric boundary extension · adaptive contexts · independent code blocks',ha='center',fontsize=11)
finish(fig,__file__,'The architecture separates transform, quantization, embedded code-block coding, context-adaptive arithmetic coding, and rate/layer organization. It distinguishes irreversible CDF 9/7 coding from reversible integer 5/3 coding. The diagram describes the stages in the chapter without implying that arbitrary byte truncation remains decodable.')
