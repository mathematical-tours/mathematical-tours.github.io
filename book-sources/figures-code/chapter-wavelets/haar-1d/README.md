# Figure 4.4

Wavelets.

## Figure 4.4

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
1-D Haar multiresolution projection $P_{V_j}f$ of a function $f$.\relax
```

Each Haar projection is the exact average on its dyadic interval. Increasing j coarsens the approximation: the number of intervals decreases from 64 to 8. The same fine signal is shown in gray throughout, and the squared projection error is checked to increase when intervals are merged.

Omitted from the current comparison PDF. Stable identifier: `wavelets--haar-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--haar-1d
```

Matching asset directory: `figures/chapter-wavelets/haar-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
