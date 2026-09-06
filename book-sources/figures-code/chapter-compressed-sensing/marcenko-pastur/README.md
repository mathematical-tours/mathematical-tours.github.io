# Figure 12.2

Compressed Sensing.

## Figure 12.2

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Empirical eigenvalue distributions approaching the Marchenko--Pastur law.\relax
```

Each histogram comes from eigenvalues of BᵀB with independent N(0,1/P) entries and the stated aspect ratio s/P. The overlay is the exact Marchenko–Pastur density for beta<1, so there is no atom at zero. Density normalization is not confused with singular-value density.

Omitted from the current comparison PDF. Stable identifier: `compressed-sensing--marcenko-pastur`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id compressed-sensing--marcenko-pastur
```

Matching asset directory: `figures/chapter-compressed-sensing/marcenko-pastur`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
