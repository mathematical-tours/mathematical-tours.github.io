# Figure 13.3

Basics of Machine Learning.

## Figure 13.3

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
PCA projections of the data into two and three dimensions.\relax
```

Both views use the same centered Iris data and orthonormal PCA axes, without standardizing the centimeter measurements. Species labels are used only for color, not to fit the projection.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--pca`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id machine-learning--pca
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/pca`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
