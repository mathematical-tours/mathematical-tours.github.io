# Figure 15.3

Optimization \& Machine Learning: Advanced Topics.

## Figure 15.3

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Objective-error evolution for batch gradient descent in logistic classification.\relax
```

Batch gradient descent is evaluated on one regularized logistic finite-sum problem. The quadratic penalty ensures a finite unique minimizer; a high-accuracy independent solve supplies a numerical objective reference. All three steps are below 2/L, and gaps below 1e-16 are displayed at the numerical floor.

Omitted from the current comparison PDF. Stable identifier: `optim-ml-advanced--bgd`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id optim-ml-advanced--bgd
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-advanced-topics/bgd`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
