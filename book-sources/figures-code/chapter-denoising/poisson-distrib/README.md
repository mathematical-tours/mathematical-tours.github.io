# Figure 7.26

Denoising.

## Figure 7.26

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.27**.

Exact current book caption (LaTeX):

```tex
Poisson distributions for various $\la $.\relax
```

Displays discrete Poisson probability masses, not continuous density curves, for means 1,4,10,20. The parameter is both the expectation and variance. The plotted support is truncated at 45 counts; omitted probability is recorded.

Omitted from the current comparison PDF. Stable identifier: `denoising--poisson-distrib`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--poisson-distrib
```

Matching asset directory: `figures/chapter-denoising/poisson-distrib`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
