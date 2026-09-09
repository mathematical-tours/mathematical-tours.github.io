from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from figure_style import *
from duality_data import chirp
from transforms import fwt,walsh_sequency
x,f=chirp();s=np.abs(fwt(f));order=np.argsort(walsh_sequency(len(f)))
fig,axs=plt.subplots(1,2,figsize=(7.1,2.65),sharey=True)
for ax,y,title,label in zip(axs,[s,s[order]],['Natural order','Sequency order'],['Character index','Sign changes']):
    ax.plot(np.arange(len(f)),y,color=TEAL,lw=.9)
    ax.set(title=title,xlim=(0,255),ylim=(0,1.07*max(s)),xticks=[0,64,128,192,255])
    axes_box(ax,label)
axs[0].set_ylabel(r'$|\mathcal{W}_8f|$')
save(fig,'chapter-duality','walsh-spectrum')
