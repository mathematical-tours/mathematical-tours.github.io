# Figure 8.6

Variational Priors and Regularization.

## Figure 8.6

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Regularized gradient norm $\sqrt { \norm {\nabla f(x)}^2 + \epsilon ^2 }$.\relax
```

All panels evaluate sqrt(|gradient f|²+epsilon²) on the same photograph and share the same 0 to .8 color range. In particular the increasing nonzero baseline remains visible instead of being removed by per-panel contrast normalization.

Omitted from the current comparison PDF. Stable identifier: `variational-priors--discrete-operators-tveps`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id variational-priors--discrete-operators-tveps
```

Matching asset directory: `figures/chapter-variational-priors-and-regularization/discrete-operators-tveps`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
