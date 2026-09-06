# Figure 4.16

Wavelets.

## Figure 4.16

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.15**.

Exact current book caption (LaTeX):

```tex
Steps of the anisotropic wavelet transform.\relax
```

Restores the flower photograph and the original row-then-column presentation. Separators follow the four-level one-dimensional packing after the first pass and its full tensor-product grid after the second pass. Every coefficient is computed, not copied from the old coefficient image. Coefficient displays use a shared signed logarithmic grayscale with zero at mid-gray; the coarse approximation is displayed at its image intensity scale. Colored boundaries follow the actual dyadic subband slices, from red at the finest level to blue, green and purple at coarser levels. Display compression does not alter the numerical coefficients.

Omitted from the current comparison PDF. Stable identifier: `wavelets--wavalgo-2d-aniso`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--wavalgo-2d-aniso
```

Matching asset directory: `figures/chapter-wavelets/wavalgo-2d-aniso`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
