"""Real and imaginary parts of four characters of Z/12Z."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'tools'))
from figure_style import *

x = np.arange(12)
t = np.linspace(0, 11, 400)
fig, axs = plt.subplots(2, 4, figsize=(8.0, 3.35), sharex=True, sharey=True)
for j in range(4):
    for row, (part, color, label) in enumerate([
        (np.cos, TEAL, 'Re'), (np.sin, BLUE, 'Im')]):
        ax = axs[row, j]
        ax.plot(t, part(2*np.pi*j*t/12), color=color, alpha=.55, lw=1)
        ax.plot(x, part(2*np.pi*j*x/12), 'o', ms=3.5, color=color)
        ax.axhline(0, color=MUTED, lw=.45, zorder=0)
        ax.set(xlim=(-.3,11.3), ylim=(-1.15,1.15), xticks=[0,5,10],
               yticks=[-1,0,1], title=rf'$\mathrm{{{label}}}\,\chi_{{{j}}}$')
        axes_box(ax, '$x$' if row else None)
characters = np.exp(2j*np.pi*np.outer(x,x)/12)
assert np.allclose(characters @ characters.conj().T, 12*np.eye(12))
save(fig, 'chapter-finite-groups', 'cyclic-characters')
