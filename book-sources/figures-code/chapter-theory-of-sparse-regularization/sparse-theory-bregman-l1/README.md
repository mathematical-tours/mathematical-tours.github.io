# Figure 11.6

Theory of Sparse Regularization.

## Figure 11.6

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Bregman divergence controls the $\ell ^1$ error on coordinates where $\eta $ does not saturate.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Explicitly sets x0=0, as required on nonsaturated coordinates for J=abs(x). For abs(eta)<1 the supporting line eta*x lies strictly below abs(x) away from zero. The lower bound is D_eta(x|0)>= (1-abs(eta))*abs(x), not a norm bound on saturated coordinates.

Omitted from the current comparison PDF. Stable identifier: `sparse-theory--sparse-theory-bregman-l1`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-theory--sparse-theory-bregman-l1
```

Matching asset directory: `figures/chapter-theory-of-sparse-regularization/sparse-theory-bregman-l1`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
