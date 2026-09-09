from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from tikz_figure import build
build(Path(__file__).with_name('figure.tex'),'chapter-duality','geometric-reciprocity')
