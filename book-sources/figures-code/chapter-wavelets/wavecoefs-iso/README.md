# Figure 4.17

Wavelets.

## Figure 4.17

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.16**.

Exact current book caption (LaTeX):

```tex
Anisotropic (left) versus isotropic (right) wavelet coefficients.\relax
```

Restores the original two packed coefficient layouts using the supplied flower image. Full horizontal and vertical separators identify independent tensor-product scales on the left; nested quadrant boundaries identify recursive low-low subdivision on the right. The filter family and display map are identical. Coefficient displays use a shared signed logarithmic grayscale with zero at mid-gray; the coarse approximation is displayed at its image intensity scale. Colored boundaries follow the actual dyadic subband slices, from red at the finest level to blue, green and purple at coarser levels. Display compression does not alter the numerical coefficients.

Omitted from the current comparison PDF. Stable identifier: `wavelets--wavecoefs-iso`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--wavecoefs-iso
```

Matching asset directory: `figures/chapter-wavelets/wavecoefs-iso`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
