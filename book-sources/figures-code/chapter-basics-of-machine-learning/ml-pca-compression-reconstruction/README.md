# Figure 13.5

Basics of Machine Learning.

## Figure 13.5

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Linear compression followed by reconstruction.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Both R and S have p rows and k columns. Compression applies S transpose; reconstruction applies R, so the product is a p-by-p map of rank at most k. Orthogonality is not assumed for arbitrary R and S; it is justified later when reducing to orthogonal projections.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--ml-pca-compression-reconstruction`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id machine-learning--ml-pca-compression-reconstruction
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/ml-pca-compression-reconstruction`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
