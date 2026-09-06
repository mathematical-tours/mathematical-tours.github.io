# Figure 13.6

Basics of Machine Learning.

## Figure 13.6

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Panel 1: the linear-program bound in~\eqref {eq-variational-pca}. Panel 2: matrix $B$. Panel 3: its orthogonal extension $\tilde B$.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The left feasible region is a triangle for p=2, k=1. The right is the unit cube truncated by sum beta <= 2 for p=3, k=2. Its cutting face is shown with nonzero visible area; dashed edges distinguish the hidden cube edges from the removed corner (1,1,1). Red dots mark valid maximizers (1,0) and (1,1,0) for ordered nonnegative eigenvalues, without asserting uniqueness. The matrix B has p rows and k orthonormal columns. The vector b_i is in R^k and appears transposed as a row. Its squared row norms sum to k, the squared Frobenius norm of B. These dimensions follow the current proof rather than ambiguous letters in the scan. The first k columns are B; the remaining p-k columns form B perpendicular. Every full row of the square orthogonal matrix has norm one. Restricting a row to the first k entries cannot increase its norm. When k=p, the complementary block is empty.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--pca-var-proof`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id machine-learning--pca-var-proof
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/pca-var-proof`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
