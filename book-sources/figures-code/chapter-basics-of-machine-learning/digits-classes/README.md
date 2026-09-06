# Figure 13.19

Basics of Machine Learning.

## Figure 13.19

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Digit classification. Left: the first nine class probabilities on an affine two-dimensional PCA slice through the training mean. Right: the mixture of all ten class colors, with projected observations overlaid.\relax
```

A ten-class multinomial logistic model is trained on 75% of the real digit images. The nine small panels show classes 0–8 on one probability scale; the RGB mixture includes all ten probabilities. Grid points are lifted from the two-PC plane back to 64 pixel coordinates before prediction. Dots are held-out images projected onto these axes and colored by their true class. Their full-space predictions need not coincide with the affine slice at the projected locations. PCA is fit on training data only.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--digits-classes`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--digits-classes
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/digits-classes`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
