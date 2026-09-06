# Figure 13.17

Basics of Machine Learning.

## Figure 13.17

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Effect of class separation on predicted classification probabilities.\relax
```

Three logistic models are fitted to paired Gaussian clouds with a shared random realization and increasing class-mean separation. Probability maps use the same 0–1 color scale and finite quadratic regularization, so complete separation does not send an unregularized coefficient to infinity.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--separation-influ`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--separation-influ
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/separation-influ`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
