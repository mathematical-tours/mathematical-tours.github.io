# Figure 13.21

Basics of Machine Learning.

## Figure 13.21

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Nonlinear classification using a Gaussian kernel.\relax
```

A Gaussian-kernel logistic classifier is fitted independently for each bandwidth on one fixed two-moons sample. Decision regions are the computed sign of its score, not hand-drawn curves or uncalibrated probabilities. The sample is synthetic and labeled as such.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--classes-kernel`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id machine-learning--classes-kernel
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/classes-kernel`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
