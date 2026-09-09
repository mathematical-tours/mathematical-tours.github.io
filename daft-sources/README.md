# The Discrete Algebra of the Fourier Transform

Gabriel Peyré — CNRS & DMA, École Normale Supérieure, PSL University.

This is an **adapted English translation** of the author's French book
[*L'algèbre discrète de la transformée de Fourier*](https://www.editions-ellipses.fr/fr/accueil/12783-23434-l-algebre-discrete-de-la-transformee-de-fourier-niveau-m1-9782729818678.html),
published by Ellipses in **2004**, ISBN **9782729818678**.
The September 2026 edition revises the mathematics and exposition and uses Python
throughout the computational examples.

- [Read the book](DAFT-EN.pdf).
- [Explore all 68 figures, their code, and Python notebooks](https://mathematical-tours.github.io/daft/figures/).
- [Major corrections](corrections.md).
- [Download the Python examples and checks](downloads/daft-python.zip).
- [Previous edition and its MATLAB/Maple sources](legacy/README.md).

## Repository layout

| Directory | Contents |
| --- | --- |
| `chapters-tex/` | Eight chapters with exercises and chapter-end solutions, two Python appendices and shared typography |
| `listings/` | Executable Python modules; `tex/` selects the exact functions printed in the book |
| `figures/chapter-title/figure-name/` | The 68 generated figure PDFs used by the book |
| `figures-code/chapter-title/figure-name/` | Each figure's Python generator and, where applicable, TikZ source |
| `data/` | The supplied flower photograph and historical cover artwork |
| `notebooks/chapter-title/` | Runnable notebooks for all 60 Python figures |
| `tests/` | Independent mathematical and numerical checks |
| `tools/` | Reproduction, verification and packaging utilities |
| `docs/` | Figure/program inventories and verification reports |
| `legacy/` | Unchanged previous edition and its SHA-256 inventory |

## Run the Python examples

Use Python 3.12. Create a virtual environment and install the recorded dependencies:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
make check PYTHON=.venv/bin/python
```

The checks exercise every function printed in all 88 numbered programs. They
compare transforms with direct sums, test algebraic identities exactly, solve
manufactured differential equations, and check coding algorithms including
nonbinary errors. `docs/program-verification.json` records the results and the
execution coverage of printed functions; it does not claim exhaustive branch coverage.

For an individual example, run Python from this directory:

```python
import sys
sys.path.insert(0, "listings")
import numpy as np
from fft_algorithms import fft_dit
from error_correcting_codes import bch_example

print(fft_dit([1, 2, 0, -1]))
code, word, received, corrected = bch_example()
```

Modules contain their own imports and shared helper functions. The book prints
compact function excerpts rather than repeating these imports in every program.
All new numerical illustrations use synthetic data or the supplied flower image.
The cover features a detail of a still-life painting by [Louis Peyré (1923–2012)](https://www.louispeyre.com/). `docs/cover-artwork.json` records the source image and checksum.

## Build the documents

A TeX Live installation with `latexmk`, pdfLaTeX, New PX, Source Sans Pro,
TikZ and the usual AMS packages is required. TeX Live 2022 was used for this edition.

```sh
make book PYTHON=.venv/bin/python
make documents PYTHON=.venv/bin/python
```

The figure PDFs are supplied, so rebuilding the book does not require regenerating
them. To reproduce all figures, or one selected figure:

```sh
make figures PYTHON=.venv/bin/python
.venv/bin/python tools/regenerate_figures.py --figure fig-cyclic-characters
```

Then rebuild the book and run `make gallery PYTHON=.venv/bin/python` to refresh
the gallery previews, figure numbers, PDF page links, notebooks, and source bundle.
The gallery reads the current captions and the compiled PDF destinations.
Rebuild the book before rebuilding the gallery whenever pagination changes.
`make python-download` rebuilds the small standalone code archive.
`make clean` removes TeX working files while retaining the deliverable PDFs.

The figure gallery uses the same vector PDFs and generation code as the book.
Its notebooks display the full plotting code and share a versioned source archive,
so downloaded notebooks continue to use their corresponding figure inputs.
