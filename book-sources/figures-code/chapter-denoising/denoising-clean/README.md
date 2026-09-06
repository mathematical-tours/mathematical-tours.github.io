# Figure 7.5

Denoising.

## Figure 7.5

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Left: clean image, center: noisy image, right: denoised image.\relax
```

The estimate is computed by hard thresholding stationary Daubechies-4 wavelet details at 2.7 sigma and applying the exact inverse stationary transform. The coarsest approximation is retained. Measured SNR is shown for the noisy and reconstructed images.

Omitted from the current comparison PDF. Stable identifier: `denoising--denoising-clean`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--denoising-clean
```

Matching asset directory: `figures/chapter-denoising/denoising-clean`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
