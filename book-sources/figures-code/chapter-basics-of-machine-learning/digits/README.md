# Figure 13.18

Basics of Machine Learning.

## Figure 13.18

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Digit images and their two- and three-dimensional PCA projections.\relax
```

Uses the real 1797-image optical handwritten-digits dataset bundled with scikit-learn. Pixel intensities are divided by 16, and both projections use one centered PCA fit. Labels color the plots but do not determine the principal components.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--digits`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id machine-learning--digits
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/digits`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
