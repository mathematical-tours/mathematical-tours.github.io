# Figure 5.14

Linear and Nonlinear Approximation.

## Figure 5.14

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Linear (top row) and nonlinear (bottom row) Fourier approximation.\relax
```

The top row keeps a fixed disk of low spatial frequencies in a real Fourier basis. The bottom row chooses the largest M real coefficients of the same transform. Each column has an identical coefficient budget, image, and display range.

Omitted from the current comparison PDF. Stable identifier: `approximation--fourier-approx-2d-linnonlin`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--fourier-approx-2d-linnonlin
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/fourier-approx-2d-linnonlin`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
