# Figure 4.2

Wavelets.

## Figure 4.2

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Haar scaling and wavelet functions. Both have amplitude $2^{-j/2}$ and unit $L^2$ norm.\relax
```

Plots the exact Haar functions at j=0,n=0. Half-open interval conventions are marked with filled and open endpoints, without drawing spurious vertical graph segments at jumps. General atoms follow by translation and multiplication by 2^(-j/2) at scale 2^j; both displayed functions have unit L2 norm.

Omitted from the current comparison PDF. Stable identifier: `wavelets--haard-scaling-wav`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--haard-scaling-wav
```

Matching asset directory: `figures/chapter-wavelets/haar-scaling-wavelet`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
