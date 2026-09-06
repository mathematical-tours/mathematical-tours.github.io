# Figure 7.30

Denoising.

## Figure 7.30

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.31**.

Exact current book caption (LaTeX):

```tex
Multiplicative Gamma noise for increasing numbers $K$ of averaged looks, with standard deviation $\sigma =K^{-1/2}$.\relax
```

Multipliers have the explicit Gamma(shape=K,scale=1/K) law, hence mean one and variance 1/K. The Flower intensity is kept strictly positive for compatibility with the subsequent log transform. All panels use one display range; saturation is a display effect rather than truncation of the simulated model.

Omitted from the current comparison PDF. Stable identifier: `denoising--multnoise`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--multnoise
```

Matching asset directory: `figures/chapter-denoising/multnoise`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
