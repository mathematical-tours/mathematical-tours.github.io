# Figure 5.16

Linear and Nonlinear Approximation.

## Figure 5.16

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Segmentation of the wavelet coefficients into regular and singular parts.\relax
```

Every dot represents an actual Haar detail coefficient. Color is determined by whether its dyadic support intersects a known jump, and area by magnitude. Periodicity adds the boundary jump. Haar is used so support intersection is exact and visible.

Omitted from the current comparison PDF. Stable identifier: `approximation--singular-support-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--singular-support-1d
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/singular-support-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
