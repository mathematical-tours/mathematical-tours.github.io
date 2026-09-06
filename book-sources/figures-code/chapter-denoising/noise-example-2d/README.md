# Figure 7.4

Denoising.

## Figure 7.4

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Two-dimensional additive noise: the clean image $f_0$ and the observation $f=f_0+w$, with $w\sim \mathcal N(0,\sigma ^2 I)$.\relax
```

A single fixed white-Gaussian realization is added to the supplied flower luminance. Both images use the range [0,1]; clipping is for display only. This same clean/noisy pair is reused in the subsequent thresholding comparisons.

Omitted from the current comparison PDF. Stable identifier: `denoising--noise-example-2d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--noise-example-2d
```

Matching asset directory: `figures/chapter-denoising/noise-example-2d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
