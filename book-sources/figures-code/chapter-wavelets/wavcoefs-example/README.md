# Figure 4.21

Wavelets.

## Figure 4.21

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.20**.

Exact current book caption (LaTeX):

```tex
Images (top) and their wavelet coefficients (bottom).\relax
```

Retains the original image-over-coefficients arrangement. A synthetic square, the supplied flower image and a central flower detail compare simple edges with photographic structure. Every lower panel is computed from the image directly above, with the same transform and display scale. Coefficient displays use a shared signed logarithmic grayscale with zero at mid-gray; the coarse approximation is displayed at its image intensity scale. Colored boundaries follow the actual dyadic subband slices, from red at the finest level to blue, green and purple at coarser levels. Display compression does not alter the numerical coefficients.

Omitted from the current comparison PDF. Stable identifier: `wavelets--wavcoefs-example`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--wavcoefs-example
```

Matching asset directory: `figures/chapter-wavelets/wavcoefs-example`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
