from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
t,f=signal_1d(1024);fig,axs=canvas(2,3,height=4.8,width=11)
for col,(fine,coarse) in enumerate([(64,32),(32,16),(16,8)]):
    p=np.repeat(f.reshape(fine,-1).mean(axis=1),len(f)//fine);q=np.repeat(f.reshape(coarse,-1).mean(axis=1),len(f)//coarse);d=p-q
    axs[0,col].plot(t,p,color=MUTED,lw=.9,label='Finer projection');axs[0,col].step(t,q,where='post',label='Coarser projection');axs[0,col].set(title=f'{fine} → {coarse} intervals',ylim=(0,1.2),xlabel=r'$t$');axs[0,col].legend(fontsize=8)
    axs[1,col].step(t,d,where='post',color=RED);axs[1,col].set(title='Detail = finer − coarser',xlabel=r'$t$',ylim=(-.35,.35))
    assert abs(np.dot(q,d))<1e-10 and np.max(abs(p-q-d))<1e-12
finish(fig,__file__,'Three successive Haar resolution pairs are arranged above their corresponding detail functions, including the additional 64-to-32 interval pair. Each column satisfies finer = coarser + detail and exact coarse/detail orthogonality. The resolution and detail rows use common vertical ranges.',checks={'reconstruction_identity':True,'coarse_detail_orthogonality':True})
