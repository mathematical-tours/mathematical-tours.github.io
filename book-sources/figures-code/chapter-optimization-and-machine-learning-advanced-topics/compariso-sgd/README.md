# Figure 15.7

Optimization \& Machine Learning: Advanced Topics.

## Figure 15.7

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Objective gaps for SGD, averaged SGD (SGA), and SAG, displayed on logarithmic axes.\relax
```

SGD, Cesàro-averaged SGD (SGA), and SAG use the same logistic problem and sampled-index stream. SGA averages its own auxiliary iterates. SAG maintains the exact mean of a table of past gradients, starts with zero entries, and updates the sampled entry before stepping; its direction is not called unbiased. Schedules and a conservative SAG step 1/(16 L_i) are explicit in the shared source.

Omitted from the current comparison PDF. Stable identifier: `optim-ml-advanced--compariso-sgd`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id optim-ml-advanced--compariso-sgd
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-advanced-topics/compariso-sgd`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
