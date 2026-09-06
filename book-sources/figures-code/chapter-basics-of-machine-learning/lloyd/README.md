# Figure 13.9

Basics of Machine Learning.

## Figure 13.9

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Continuous $k$-means (Lloyd) iterations $0,2,3,5,30$. The top row uses a uniform density; the bottom row uses the nonuniform density shown in grayscale.\relax
```

Continuous Lloyd updates are approximated by weighted quadrature on a 90 by 90 grid. Both rows start from the same 28 centers. Nonuniform centroids use density weights, not the unweighted mean of the Voronoi pixels. Cell boundaries are the exact Euclidean Voronoi edges of the current centroids.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--lloyd`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--lloyd
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/lloyd`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
