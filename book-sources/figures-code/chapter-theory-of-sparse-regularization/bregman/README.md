# Figure 11.5

Theory of Sparse Regularization.

## Figure 11.5

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Visualization of Bregman divergences.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The left panel uses a strictly convex quadratic and its tangent; the positive vertical gap is the Bregman divergence. The right panel uses J=abs(x), with x0 and x on the same positive affine branch, so the divergence vanishes even when x differs from x0.

Omitted from the current comparison PDF. Stable identifier: `sparse-theory--bregman`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-theory--bregman
```

Matching asset directory: `figures/chapter-theory-of-sparse-regularization/bregman`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
