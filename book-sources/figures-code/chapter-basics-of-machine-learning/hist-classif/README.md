# Figure 13.14

Basics of Machine Learning.

## Figure 13.14

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Classification boundaries produced by $R$-nearest neighbors.\relax
```

Nearest-neighbor classifiers are fitted in the displayed two-dimensional Iris PCA coordinates. Exact majority votes on a fine grid generate the three decision maps. These are illustrative training-data boundaries, not a held-out accuracy estimate.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--hist-classif`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--hist-classif
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/hist-classif`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
