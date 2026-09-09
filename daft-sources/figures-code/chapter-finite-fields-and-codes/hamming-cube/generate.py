from pathlib import Path
import sys
root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(root/'tools'))
from tikz_figure import build
build(Path(__file__).with_name('figure.tex'), 'chapter-finite-fields-and-codes', 'hamming-cube')
