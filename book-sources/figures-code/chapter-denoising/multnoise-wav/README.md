# Figure 7.32

Denoising.

## Figure 7.32

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.33**.

Exact current book caption (LaTeX):

```tex
Left: noisy image, center: denoising after variance stabilization, right: denoising without variance stabilization.\relax
```

The stabilized target is log(f0): the known mean log multiplier is subtracted before thresholding, and no extra centering constant is added during exponentiation. Both methods use identical noisy data and independently oracle-tuned thresholds. Exponentiating the estimate may introduce bias; the displayed SNR values are actual results, not a claim of unbiasedness or universal superiority.

Omitted from the current comparison PDF. Stable identifier: `denoising--multnoise-wav`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--multnoise-wav
```

Matching asset directory: `figures/chapter-denoising/multnoise-wav`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
