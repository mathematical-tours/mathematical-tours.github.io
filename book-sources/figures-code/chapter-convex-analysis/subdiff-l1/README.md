# Figure 18.3

Convex Analysis.

## Figure 18.3

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Subdifferential of the absolute value and a piecewise affine convex function.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. At each kink the entire closed interval of intermediate slopes is drawn, not just its endpoints. The second function uses representative slopes m1<m2<0<m3 at knots a<0 and 0; the original does not specify numerical slopes. The lower panels are set-valued graphs, so their vertical segments are essential.

Omitted from the current comparison PDF. Stable identifier: `convex-analysis--subdiff-l1`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id convex-analysis--subdiff-l1
```

Matching asset directory: `figures/chapter-convex-analysis/subdiff-l1`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
