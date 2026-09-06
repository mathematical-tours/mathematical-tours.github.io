# Figure 8.9

Variational Priors and Regularization.

## Figure 8.9

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
SNR as a function of time $t$ for flows (top) and $\la $ for regularization (bottom).\relax
```

The curves and marked grid maxima are measured on the exact reconstructions used by Figure 8.8. Each axis states whether the parameter is time or regularization weight, and logarithmic abscissae resolve small values. All methods share one noise realization.

Omitted from the current comparison PDF. Stable identifier: `variational-priors--denoising-cameraman-snr`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id variational-priors--denoising-cameraman-snr
```

Matching asset directory: `figures/chapter-variational-priors-and-regularization/denoising-cameraman-snr`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
