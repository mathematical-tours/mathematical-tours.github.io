# Figure 8.1

Variational Priors and Regularization.

## Figure 8.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Discrete gradient vectors.\relax
```

Vectors use the actual forward finite differences, with columns drawn horizontally and rows vertically in image coordinates. The arrows point toward increasing intensity; the two detail panels recompute no new field and retain the same physical arrow scale.

Omitted from the current comparison PDF. Stable identifier: `variational-priors--grad-discrete`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id variational-priors--grad-discrete
```

Matching asset directory: `figures/chapter-variational-priors-and-regularization/grad-discrete`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
