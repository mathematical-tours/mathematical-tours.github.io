# Figure 5.15

Linear and Nonlinear Approximation.

## Figure 5.15

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Singular and regular regions of a signal (left) and an image (right).\relax
```

Red marks identify neighborhoods of singular points in one dimension and a discontinuity curve in two dimensions. The underlying functions are piecewise smooth; the shaded jump neighborhoods represent supports that intersect singularities, not an additional discontinuity of the function.

Omitted from the current comparison PDF. Stable identifier: `approximation--singular-part`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--singular-part
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/singular-part`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
