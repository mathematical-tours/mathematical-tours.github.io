# Figure 7.9

Denoising.

## Figure 7.9

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Gaussian filtering at the SNR-optimal widths. From left to right: noisy signal, filtered signal, noisy image, and filtered image.\relax
```

Combines the previously separate signal and image comparisons on one row. Both oracle filter widths use exactly the observations and width grid of the preceding SNR plot. Clean references are used for selection and SNR measurement; the image panels share one fixed intensity range.

Omitted from the current comparison PDF. Stable identifier: `denoising--filtering-optimal-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--filtering-optimal-1d
```

Matching asset directory: `figures/chapter-denoising/filtering-optimal-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
