# Figure 10.2

Sparse Regularization.

## Figure 10.2

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
$\ell ^q$ balls $\enscond {x}{J_q(x) \leq 1}$ for varying $q$.\relax
```

Exact radial parameterizations of the two-dimensional ell-q unit sets preserve equal axis scales. The q=.5 set is nonconvex; q=1 is a diamond, q=1.5 and q=4 have intermediate strictly convex boundaries, q=2 is a disk, and q=infinity is a square.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--sparsity-lq`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id sparse-regularization--sparsity-lq
```

Matching asset directory: `figures/chapter-sparse-regularization/sparsity-lq`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
