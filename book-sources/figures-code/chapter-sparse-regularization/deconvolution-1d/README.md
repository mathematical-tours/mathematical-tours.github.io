# Figure 10.7

Sparse Regularization.

## Figure 10.7

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
One-dimensional deconvolution with sparsity in an orthogonal wavelet basis.\relax
```

One-dimensional deconvolution solves .5||H Psi c-y||²+.005||c||1 in a periodic orthonormal db2 basis. The synthesis coefficients, not image values, are soft-thresholded by the optimizer. The same filter produces the observation and the forward/adjoint operators.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--deconvolution-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id sparse-regularization--deconvolution-1d
```

Matching asset directory: `figures/chapter-sparse-regularization/deconvolution-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
