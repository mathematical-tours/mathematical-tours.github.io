# Figure 7.22

Denoising.

## Figure 7.22

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.23**.

Exact current book caption (LaTeX):

```tex
SNR as a function of $T/\si $ for Stein thresholding.\relax
```

Stationary Daubechies-4 thresholding is evaluated with Stein, hard, and soft nonlinearities on identical data. SNR values and grid maxima are recomputed, so the figure records the result of this experiment rather than asserting that one method is universally best.

Omitted from the current comparison PDF. Stable identifier: `denoising--stein-optimal-curve`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--stein-optimal-curve
```

Matching asset directory: `figures/chapter-denoising/stein-optimal-curve`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
