# Figure 13.8

Basics of Machine Learning.

## Figure 13.8

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Left: iterations of the $k$-means algorithm. Right: class histograms after the $k$-means optimization.\relax
```

Lloyd iterations run in the original four-dimensional feature space; PCA is only the display projection. Right-hand bars count the true species inside the final fitted clusters, preserving the distinction between unsupervised cluster IDs and class labels.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--kmeans`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id machine-learning--kmeans
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/kmeans`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
