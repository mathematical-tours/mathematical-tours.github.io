# Figure 15.1

Optimization \& Machine Learning: Advanced Topics.

## Figure 15.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
$\ell ^q$ balls $\enscond {x}{\sum _k |x_k|^q \leq 1}$ for varying $q$.\relax
```

Exact radial parameterizations of the two-dimensional ell-q unit sets preserve equal axis scales. The q=.5 set is nonconvex; q=1 is a diamond, q=1.5 and q=4 have intermediate strictly convex boundaries, q=2 is a disk, and q=infinity is a square.

Omitted from the current comparison PDF. Stable identifier: `optim-ml-advanced--ml-sparsity-lq`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id optim-ml-advanced--ml-sparsity-lq
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-advanced-topics/ml-sparsity-lq`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
