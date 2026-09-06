# Figure 10.8

Sparse Regularization.

## Figure 10.8

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Image deconvolution.\relax
```

Squared norm, Sobolev and orthogonal wavelet reconstructions use the exact grid plotted in Figure 10.9. The last panel runs coefficient-space synthesis FISTA with a normalized stationary Parseval frame; it does not pretend that analysis-threshold-synthesis is the proximal map of an analysis penalty. Each stationary coefficient is penalized by its atom L2 norm, equivalently using a unit-norm synthesis dictionary. Its lambda equals the orthogonal oracle choice and is stated as such. The frame adjoint identity is checked.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--deconvolution-2d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-regularization--deconvolution-2d
```

Matching asset directory: `figures/chapter-sparse-regularization/deconvolution-2d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
