# Figure 5.17

Linear and Nonlinear Approximation.

## Figure 5.17

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
1-D wavelet approximation.\relax
```

Nonlinear orthonormal wavelet approximations retain the largest 16,32,64,128 coefficients of the same 1024-sample signal. Fine details survive primarily near jumps.

Omitted from the current comparison PDF. Stable identifier: `approximation--wav-approx-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--wav-approx-1d
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/wav-approx-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
