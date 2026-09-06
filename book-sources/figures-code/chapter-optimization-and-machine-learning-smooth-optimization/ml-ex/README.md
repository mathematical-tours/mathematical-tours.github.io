# Figure 14.1

Optimization \& Machine Learning: Smooth Optimization.

## Figure 14.1

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Panels 1 and 2: linear regression and a linear classifier. Panel 3: the 0-1 loss and its normalized logistic and hinge upper bounds, both tight at zero margin.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The fifteen training pairs form an illustrative scalar-feature dataset. The displayed slope is their exact least-squares minimizer, x-star=sum(a_i y_i)/sum(a_i squared), approximately 0.4589392. The red vertical segment is the signed residual y_i−〈a_i,x-star〉 for the marked observation; no intercept is fitted. Positive labels lie on the positive-score side and negative labels on the negative-score side. The normal vector is perpendicular to the separating line, with an explicit right-angle marker. The feature-space origin is marked because this homogeneous decision hyperplane passes through it. Divided the logistic loss by log(2), giving ell(u)=log(1+exp(u))/log(2). This is an upper bound on the displayed zero-one loss and meets it at u=0, where both equal one. Zero margin counts as an error. The hinge loss also equals one there. The chapter already uses this normalized loss; its derivative is the sigmoid divided by log(2).

Omitted from the current comparison PDF. Stable identifier: `optim-ml-smooth--ml-ex`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id optim-ml-smooth--ml-ex
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-smooth-optimization/ml-ex`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
