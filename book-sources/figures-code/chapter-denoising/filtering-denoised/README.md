# Figure 7.6

Denoising.

## Figure 7.6

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Gaussian filtering: noisy signal, kernel, and filtered signal (top), with their Fourier magnitudes below.\relax
```

The top row contains signal, kernel and filtered signal; the bottom row contains their corresponding Fourier spectra. All six panels belong to one cyclic convolution g=y*h, G=YH. The Gaussian has unit discrete mass; the plotting floor changes only logarithmic display.

Omitted from the current comparison PDF. Stable identifier: `denoising--filtering-denoised`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--filtering-denoised
```

Matching asset directory: `figures/chapter-denoising/filtering-denoised`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
