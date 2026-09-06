# Figure 8.4

Variational Priors and Regularization.

## Figure 8.4

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Discrete Laplacian (left) and negative discrete TV gradient (right).\relax
```

Both operators use backward divergence of forward differences. The second panel explicitly states epsilon=.01: the unsmoothed quotient would be undefined at zero gradients. Each panel has a labeled signed color scale; differing amplitudes are not concealed.

Omitted from the current comparison PDF. Stable identifier: `variational-priors--discrete-operators-lapl`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id variational-priors--discrete-operators-lapl
```

Matching asset directory: `figures/chapter-variational-priors-and-regularization/discrete-operators-lapl`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
