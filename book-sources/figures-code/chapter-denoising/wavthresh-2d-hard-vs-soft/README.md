# Figure 7.12

Denoising.

## Figure 7.12

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.13**.

Exact current book caption (LaTeX):

```tex
SNR as a function of $T/\si $ for hard and soft thresholding.\relax
```

Hard and soft thresholding use the same noisy flower, orthonormal Daubechies-4 basis, and threshold grid. Each SNR value is measured from an actual reconstructed image. The maxima are reference-dependent experimental optima, not general rankings of the shrinkage rules.

Omitted from the current comparison PDF. Stable identifier: `denoising--wavthresh-2d-hard-vs-soft`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--wavthresh-2d-hard-vs-soft
```

Matching asset directory: `figures/chapter-denoising/wavthresh-2d-hard-vs-soft`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
