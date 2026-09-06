# Figure 3.4

Shannon Coding Theory.

## Figure 3.4

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Entropy is minimized at a point mass and maximized at the uniform distribution. Left: level sets of the entropy formula in the positive quadrant and the constraint $p_1+p_2=1$. Middle: its restriction to that constraint. Right: entropy level sets on the three-symbol probability simplex.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The new left panel shows actual level sets of -p1 log2(p1)-p2 log2(p2) in the positive quadrant, with 0 log2(0)=0. This formula extends beyond probability vectors; the red segment p1+p2=1 is the probability constraint. Its constrained maximum is at (1/2,1/2), where H=1 and the contour is tangent to the segment. The unconstrained maximum of the extension is instead at (1/e,1/e), explaining why the contours are centered off the constraint. The middle panel restricts entropy to the segment. The right panel retains the three-symbol simplex, where the uniform entropy is log2(3).

Omitted from the current comparison PDF. Stable identifier: `shannon-coding--entropy-extrema`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id shannon-coding--entropy-extrema
```

Matching asset directory: `figures/chapter-shannon-coding-theory/entropy-extrema`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
