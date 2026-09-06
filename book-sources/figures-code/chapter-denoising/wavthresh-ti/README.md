# Figure 7.17

Denoising.

## Figure 7.17

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.18**.

Exact current book caption (LaTeX):

```tex
SNR as a function of $T/\si $ for translation-invariant thresholding.\relax
```

Every point is measured after thresholding the stationary Daubechies-4 coefficients and reconstructing the image. Hard, soft, and Stein rules share the same observation, levels, and threshold grid. The redundant transform retains all translations; no old SNR curves are reused.

Omitted from the current comparison PDF. Stable identifier: `denoising--wavthresh-ti`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--wavthresh-ti
```

Matching asset directory: `figures/chapter-denoising/wavthresh-ti`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
