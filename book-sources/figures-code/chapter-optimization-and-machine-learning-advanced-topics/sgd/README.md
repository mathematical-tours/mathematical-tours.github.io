# Figure 15.6

Optimization \& Machine Learning: Advanced Topics.

## Figure 15.6

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Objective-error evolution for SGD in logistic classification; the dashed curve shows batch gradient descent.\relax
```

One full-gradient update costs n individual gradients and is compared with n sampled SGD updates per horizontal unit. Both optimize exactly the same quadratically regularized logistic objective, with gaps measured against one high-accuracy numerical reference. SGD uses independent uniform sampling and tau=(1/L_i)/(1+k/(3n)).

Omitted from the current comparison PDF. Stable identifier: `optim-ml-advanced--sgd`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id optim-ml-advanced--sgd
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-advanced-topics/sgd`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
