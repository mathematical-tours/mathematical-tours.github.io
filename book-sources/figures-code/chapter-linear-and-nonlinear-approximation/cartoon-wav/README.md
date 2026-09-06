# Figure 5.19

Linear and Nonlinear Approximation.

## Figure 5.19

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Wavelet approximation of a cartoon image.\relax
```

A computed cartoon approximation and its actual orthogonal wavelet coefficients show concentration around edge curves. The signed reconstruction is not clipped before SNR measurement, and the coefficient display is explicitly logarithmic. Coefficient displays use a shared signed logarithmic grayscale with zero at mid-gray; the coarse approximation is displayed at its image intensity scale. Colored boundaries follow the actual dyadic subband slices, from red at the finest level to blue, green and purple at coarser levels. Display compression does not alter the numerical coefficients.

Omitted from the current comparison PDF. Stable identifier: `approximation--cartoon-wav`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--cartoon-wav
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/cartoon-wav`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
