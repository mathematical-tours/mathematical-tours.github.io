# Figure 7.10

Denoising.

## Figure 7.10

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.11**.

Exact current book caption (LaTeX):

```tex
Denoising using thresholding of wavelet coefficients.\relax
```

Restores the original image-to-coefficients-to-thresholded-coefficients-to-image sequence using the supplied flower photograph. All detail bands use exact hard thresholding, with equality mapped to zero; the coarse band is retained. Both coefficient displays use identical scaling and boundaries. Coefficient displays use a shared signed logarithmic grayscale with zero at mid-gray; the coarse approximation is displayed at its image intensity scale. Colored boundaries follow the actual dyadic subband slices, from red at the finest level to blue, green and purple at coarser levels. Display compression does not alter the numerical coefficients.

Omitted from the current comparison PDF. Stable identifier: `denoising--wavthresh`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--wavthresh
```

Matching asset directory: `figures/chapter-denoising/wavthresh`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
