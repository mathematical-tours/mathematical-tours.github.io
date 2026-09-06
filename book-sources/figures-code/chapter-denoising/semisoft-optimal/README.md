# Figure 7.21

Denoising.

## Figure 7.21

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.22**.

Exact current book caption (LaTeX):

```tex
Left: SNR as a function of $\mu $ and $T/\si $. Right: SNR versus $\mu $, with $T/\si $ optimized separately for each $\mu $.\relax
```

The heatmap is a genuine two-parameter sweep of semi-soft orthogonal wavelet denoising on one fixed noisy image. The right curve takes a separate maximum over the threshold grid at each mu; it is not a slice through a single globally selected threshold. Tuning uses the clean reference and is explicitly oracle-based.

Omitted from the current comparison PDF. Stable identifier: `denoising--semisoft-optimal`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--semisoft-optimal
```

Matching asset directory: `figures/chapter-denoising/semisoft-optimal`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
