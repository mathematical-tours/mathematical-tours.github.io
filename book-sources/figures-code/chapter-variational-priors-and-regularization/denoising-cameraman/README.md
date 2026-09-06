# Figure 8.8

Variational Priors and Regularization.

## Figure 8.8

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Denoising using PDE flows and regularization.\relax
```

Each method is computed from the same observation and uses the clean reference only for explicitly oracle-based parameter selection. The four variational grid searches are exactly those plotted in Figure 8.9. TV regularization solves the nonsmooth ROF objective by a primal-dual algorithm; the TV flow uses epsilon=.01.

Omitted from the current comparison PDF. Stable identifier: `variational-priors--denoising-cameraman`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id variational-priors--denoising-cameraman
```

Matching asset directory: `figures/chapter-variational-priors-and-regularization/denoising-cameraman`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
