# Figure 13.13

Basics of Machine Learning.

## Figure 13.13

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Nearest neighbors.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The arrows identify the four closest observations in increasing distance order. The pale disk is the selected neighborhood; its radius is exactly the distance to the fourth observation on its boundary. Every remaining point is farther away. The panel explains the distance-selection step before voting over class labels.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--ml-nearest-neighbors`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id machine-learning--ml-nearest-neighbors
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/ml-nearest-neighbors`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
