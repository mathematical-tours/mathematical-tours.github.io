# Figure 5.1

Linear and Nonlinear Approximation.

## Figure 5.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Linear versus nonlinear wavelet approximation.\relax
```

Linear and nonlinear approximations use the same orthonormal Daubechies-4 basis and exactly 1024 of 65536 real coefficients. The linear mask is the fixed 32 by 32 coarse block; nonlinear selection sorts magnitudes with a deterministic tie rule.

Omitted from the current comparison PDF. Stable identifier: `approximation--approximation`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--approximation
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/approximation`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
