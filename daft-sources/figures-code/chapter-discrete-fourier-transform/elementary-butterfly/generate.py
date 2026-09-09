from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from tikz_figure import build
build(Path(__file__).with_name('figure.tex'),'chapter-discrete-fourier-transform','elementary-butterfly')
