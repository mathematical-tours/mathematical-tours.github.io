# Figure 13.16

Basics of Machine Learning.

## Figure 13.16

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Classification losses as functions of the negative margin $s=-y\dotp {x}{\be }$, with zero margin counted as an error. The unscaled logistic loss need not upper-bound the 0--1 loss; division by $\log 2$ gives the upper-bound normalization used above.\relax
```

Losses are evaluated exactly as functions of negative margin. The binary loss equals one at zero; open/closed dots show the discontinuity. Both logistic normalizations are distinguished, and division by log(2) makes the binary-loss upper bound tight at zero. Vertical truncation is a stated viewing range.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--losses`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id machine-learning--losses
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/losses`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
