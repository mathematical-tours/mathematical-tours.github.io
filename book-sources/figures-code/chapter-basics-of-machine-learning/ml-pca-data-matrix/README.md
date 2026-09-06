# Figure 13.4

Basics of Machine Learning.

## Figure 13.4

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Rows of the centered data matrix are observations; columns are features.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The reconstruction follows the current n-observations, p-features convention. An observation x_i is a column vector in R^p; the corresponding matrix row is its transpose. Centering is explicit, and the covariance is X transpose X divided by n.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--ml-pca-data-matrix`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--ml-pca-data-matrix
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/ml-pca-data-matrix`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
