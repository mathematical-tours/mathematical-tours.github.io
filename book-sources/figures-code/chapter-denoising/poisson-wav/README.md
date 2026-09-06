# Figure 7.29

Denoising.

## Figure 7.29

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.30**.

Exact current book caption (LaTeX):

```tex
Left: noisy image, center: denoising without variance stabilization, right: denoising after variance stabilization.\relax
```

Both denoisers use the same Poisson counts, stationary wavelet transform, and oracle threshold grid. Direct denoising uses one global count-noise scale sqrt(mean(counts)); stabilized denoising uses Anscombe’s approximate unit scale. The latter projects onto the inverse domain before applying z^2/4-3/8. This direct inverse is not claimed to be unbiased.

Omitted from the current comparison PDF. Stable identifier: `denoising--poisson-wav`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--poisson-wav
```

Matching asset directory: `figures/chapter-denoising/poisson-wav`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
