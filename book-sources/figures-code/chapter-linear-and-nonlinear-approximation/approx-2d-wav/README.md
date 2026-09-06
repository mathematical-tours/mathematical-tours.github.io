# Figure 5.18

Linear and Nonlinear Approximation.

## Figure 5.18

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
2-D wavelet approximation.\relax
```

Largest-coefficient db4 reconstructions at increasing budgets use the same acquired flower photograph. SNR is measured on the unclipped reconstruction; display uses the common unit intensity range.

Omitted from the current comparison PDF. Stable identifier: `approximation--approx-2d-wav`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--approx-2d-wav
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/approx-2d-wav`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
