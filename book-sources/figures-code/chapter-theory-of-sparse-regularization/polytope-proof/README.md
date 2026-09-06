# Figure 11.2

Theory of Sparse Regularization.

## Figure 11.2

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Geometry of the polytope characterization of exact recovery.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Uses the normalized polytope AB and q=Ax0/r, r=||x0||1, consistently with the corrected proof. The radial extension (1+t)q is constructed exactly on the ray through q and on the polygon boundary. The second panel shows h as a directed displacement to q+h inside an interior ball. Both panels illustrate failure of optimality through interior membership; no uniqueness claim is implied.

Omitted from the current comparison PDF. Stable identifier: `sparse-theory--polytope-proof`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id sparse-theory--polytope-proof
```

Matching asset directory: `figures/chapter-theory-of-sparse-regularization/polytope-proof`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
