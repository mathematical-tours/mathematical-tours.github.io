# Figure 5.2

Linear and Nonlinear Approximation.

## Figure 5.2

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Ordered coefficient magnitudes and the threshold used for nonlinear approximation. \relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The plotted values are the nonincreasing discrete sequence d_m=3.2 m^(-0.65), with M=5. T lies between d_6 and d_5, so strict thresholding retains exactly five coefficients. The source diagram incorrectly put T at the retained value d_M. The tail-energy formula uses orthonormality.

Omitted from the current comparison PDF. Stable identifier: `approximation--decay-coefs`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--decay-coefs
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/decay-coefs`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
