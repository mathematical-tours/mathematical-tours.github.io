# Figure 4.11

Wavelets.

## Figure 4.11

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Pyramidal computation of wavelet coefficients: approximation vectors $a_j$ on the top row and detail vectors $d_j$ below.\relax
```

Approximation outputs a_j are aligned on the top row; each extracted d_j is directly below the approximation from the same split. Only the approximation is recursively decomposed, and every detail band is stored once.

Omitted from the current comparison PDF. Stable identifier: `wavelets--filterbank-pyramidal`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--filterbank-pyramidal
```

Matching asset directory: `figures/chapter-wavelets/filterbank-pyramidal`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
