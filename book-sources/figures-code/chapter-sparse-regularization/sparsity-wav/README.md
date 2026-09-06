# Figure 10.1

Sparse Regularization.

## Figure 10.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Most wavelet coefficients of a natural image have small magnitude.\relax
```

The coefficient image and magnitude decay are both computed from the displayed flower. Logarithmic axes reveal many small coefficients without implying that they are exactly zero. Coefficient displays use a shared signed logarithmic grayscale with zero at mid-gray; the coarse approximation is displayed at its image intensity scale. Colored boundaries follow the actual dyadic subband slices, from red at the finest level to blue, green and purple at coarser levels. Display compression does not alter the numerical coefficients.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--sparsity-wav`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-regularization--sparsity-wav
```

Matching asset directory: `figures/chapter-sparse-regularization/sparsity-wav`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
