# Figure 7.7

Denoising.

## Figure 7.7

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Gaussian filtering with increasing width: three signals (top) and five images (bottom).\relax
```

Three signal views form the first row and five image views the second. Each row reuses a single noisy observation while increasing the standard deviation of the unit-mass periodic Gaussian filter.

Omitted from the current comparison PDF. Stable identifier: `denoising--filtering-progression`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--filtering-progression
```

Matching asset directory: `figures/chapter-denoising/filtering-progression`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
