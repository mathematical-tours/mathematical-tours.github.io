# Figure 8.3

Variational Priors and Regularization.

## Figure 8.3

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Heat flow (top) and TV flow (bottom): images $f_t$ at increasing times $t$.\relax
```

Exact discrete heat evolution is compared with explicit smoothed-TV gradient flow (epsilon=.01, step=.002 < epsilon/4). Each row states its own times because the two energy scalings differ. The common periodic flower image and unit display range expose diffusion across edges.

Omitted from the current comparison PDF. Stable identifier: `variational-priors--flow`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id variational-priors--flow
```

Matching asset directory: `figures/chapter-variational-priors-and-regularization/flow`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
