# Figure 10.10

Sparse Regularization.

## Figure 10.10

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Inpainting with Sobolev and wavelet sparsity priors.\relax
```

Reconstructs one masked image with Sobolev interpolation and two coefficient-space ell-1 synthesis models. Wavelet solutions use a small positive lambda and hence have a measured nonzero data residual; they are not labeled exact interpolants. The stationary penalty weights coefficients by their atom L2 norms, equivalent to a dictionary of unit-norm atoms. The stationary frame uses normalized analysis and its verified adjoint synthesis.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--inpainting-lena`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-regularization--inpainting-lena
```

Matching asset directory: `figures/chapter-sparse-regularization/inpainting-lena`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
