# Figure 11.1

Theory of Sparse Regularization.

## Figure 11.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Lasso solution paths $\lambda \mapsto x_\lambda $ for sparsities $s=3,6,13$, with the same measurement matrix and noise realization.\relax
```

Three LASSO homotopy paths use sparsities s=3,6,13, with one column-normalized sensing matrix, nested supports, and one noise realization. The library normalization by P is undone so the horizontal lambda matches .5||Ax-y||²+lambda||x||1. Active and inactive KKT conditions are checked at every path knot.

Omitted from the current comparison PDF. Stable identifier: `sparse-theory--homotopy`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-theory--homotopy
```

Matching asset directory: `figures/chapter-theory-of-sparse-regularization/homotopy`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
