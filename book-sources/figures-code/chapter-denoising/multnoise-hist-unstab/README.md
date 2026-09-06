# Figure 7.31

Denoising.

## Figure 7.31

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.32**.

Exact current book caption (LaTeX):

```tex
Histogram of multiplicative noise before (left) and after (right) stabilization.\relax
```

The histograms use the same Gamma multipliers before and after logarithmic stabilization. Centering uses c=digamma(K)-log(K), and the exact transformed density includes the exponential Jacobian. A Gaussian with variance trigamma(K) is shown only as an approximation; log-Gamma noise is not exactly Gaussian.

Omitted from the current comparison PDF. Stable identifier: `denoising--multnoise-hist-unstab`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--multnoise-hist-unstab
```

Matching asset directory: `figures/chapter-denoising/multnoise-hist-unstab`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
