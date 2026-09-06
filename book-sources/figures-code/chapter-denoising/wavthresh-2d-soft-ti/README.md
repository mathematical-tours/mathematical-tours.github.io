# Figure 7.16

Denoising.

## Figure 7.16

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.17**.

Exact current book caption (LaTeX):

```tex
Soft thresholding in an orthogonal wavelet basis (left) and translation-invariant hard thresholding (right).\relax
```

Compares orthogonal soft thresholding with stationary hard thresholding on the same noisy flower. Both use Daubechies-4 filters, four levels, and the same oracle threshold grid. The stationary transform uses every translation and is inverted exactly; it is not an average over only a few shifts.

Omitted from the current comparison PDF. Stable identifier: `denoising--wavthresh-2d-soft-ti`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--wavthresh-2d-soft-ti
```

Matching asset directory: `figures/chapter-denoising/wavthresh-2d-soft-ti`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
