# Figure 3.2

Shannon Coding Theory.

## Figure 3.2

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
The concave function $h(u)=-u\log _2u$ underlying entropy.\relax
```

Plots h(u)=-u log_2 u with its continuous value zero at u=0. The interior maximum is marked at u=1/e and height 1/(e ln 2); the curve is strictly concave for u>0.

Omitted from the current comparison PDF. Stable identifier: `shannon-coding--the-concave-function-underlying-entropy`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id shannon-coding--the-concave-function-underlying-entropy
```

Matching asset directory: `figures/chapter-shannon-coding-theory/the-concave-function-underlying-entropy`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
