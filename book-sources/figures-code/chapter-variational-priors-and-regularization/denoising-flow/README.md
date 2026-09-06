# Figure 8.7

Variational Priors and Regularization.

## Figure 8.7

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Denoised images $f_t$ at several times $t$: Sobolev flow (top) and TV flow (bottom).\relax
```

Flow snapshots use one fixed noisy flower. Heat and smoothed-TV times are stated separately, and SNR is recomputed before display clipping. TV time integration uses a stable step below epsilon/4.

Omitted from the current comparison PDF. Stable identifier: `variational-priors--denoising-flow`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id variational-priors--denoising-flow
```

Matching asset directory: `figures/chapter-variational-priors-and-regularization/denoising-flow`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
