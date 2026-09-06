# Figure 10.9

Sparse Regularization.

## Figure 10.9

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Reconstruction SNR as a function of the regularization parameter $\lambda $, swept logarithmically from $10^{-5}$ to $1$. Red dots mark the measured maxima.\relax
```

Measured deconvolution SNR curves use the same observation, parameters, and solved objectives as Figure 10.8. Red markers are oracle maxima on the finite lambda grid, not an automatic parameter-selection method available without a clean reference.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--deconvolution-2d-snr`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-regularization--deconvolution-2d-snr
```

Matching asset directory: `figures/chapter-sparse-regularization/deconvolution-2d-snr`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
