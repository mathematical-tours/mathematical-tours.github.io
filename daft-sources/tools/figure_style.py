"""Shared plotting conventions for the DAFT figures."""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'build' / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TEAL, RED, BLUE, INK, MUTED = '#146c78', '#c34e42', '#3478ab', '#203746', '#667780'
plt.rcParams.update({
    'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 10,
    'axes.titlesize': 11, 'axes.labelsize': 10, 'xtick.labelsize': 9,
    'ytick.labelsize': 9, 'axes.linewidth': .65, 'axes.edgecolor': MUTED,
    'axes.labelcolor': INK, 'text.color': INK, 'xtick.color': MUTED,
    'ytick.color': MUTED, 'lines.linewidth': 1.5, 'pdf.fonttype': 42,
    'ps.fonttype': 42, 'figure.constrained_layout.use': True,
    'savefig.facecolor': 'white', 'figure.facecolor': 'white',
})

def axes_box(ax, xlabel=None, ylabel=None):
    """Keep all four spines, as in the book's boxed-axis convention."""
    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.tick_params(direction='in', top=True, right=True, length=3, width=.6)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    return ax

def discrete(ax, x, y, color=TEAL):
    marker, stems, baseline = ax.stem(x, y, basefmt=' ')
    marker.set(color=color, marker='o', markersize=3.3)
    stems.set(color=color, linewidth=.8)
    return marker

def save(fig, chapter, name):
    folder = ROOT / 'figures' / chapter / name
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / 'figure.pdf', metadata={'Creator': 'DAFT Python figures',
                'CreationDate': None, 'ModDate': None}, bbox_inches='tight')
    preview = ROOT / 'build' / 'figure-previews' / chapter
    preview.mkdir(parents=True, exist_ok=True)
    fig.savefig(preview / f'{name}.png', dpi=170, bbox_inches='tight')
    plt.close(fig)
