# Figure 13.7

Basics of Machine Learning.

## Figure 13.7

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
$k$-means clusters according to Voronoi cells.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Cell edges are the exact perpendicular bisectors of the displayed centroids. The dashed centroid segment and right-angle marker make that construction explicit. Colored observations lie in their nearest-centroid cells; the common vertex is equidistant from all three centroids. This illustrates assignment, not a claim that these centroids are already the cluster means. Boundary ties require the convention in the text.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--ml-voronoi`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--ml-voronoi
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/ml-voronoi`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
