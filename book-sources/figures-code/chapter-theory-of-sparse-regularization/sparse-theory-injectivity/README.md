# Figure 11.4

Theory of Sparse Regularization.

## Figure 11.4

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Coefficient trajectories along a support-reducing direction.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Three affine coefficient trajectories have a marked interval with both endpoints shown. The first coefficient vanishes at the left endpoint and the third at the right. The figure states that x is a minimizer: this assumption forces the affine l1 norm to have zero slope while signs stay fixed. Constant fidelity follows from Ah=0. The norm remains constant at the endpoints by continuity, though the sign pattern changes there.

Omitted from the current comparison PDF. Stable identifier: `sparse-theory--sparse-theory-injectivity`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id sparse-theory--sparse-theory-injectivity
```

Matching asset directory: `figures/chapter-theory-of-sparse-regularization/sparse-theory-injectivity`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
