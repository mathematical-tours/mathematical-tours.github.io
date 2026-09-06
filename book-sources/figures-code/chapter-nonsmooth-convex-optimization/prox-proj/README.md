# Figure 19.1

Nonsmooth Convex Optimization.

## Figure 19.1

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Proximal map and projection map.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The left panel uses a quarter disk; all displayed input-to-projection segments lie in the appropriate outward normal cones. The right panel explicitly names f(u,v)=(u^2+2v^2)/2 and uses its nested sublevel sets with the exact proximal pair x=(1.4,1.35), z=(0.7,0.45). The highlighted ellipse radii are computed from f(z), so z lies exactly on its boundary. The relevant sublevel set depends on the input; a fixed set does not represent the entire proximal map.

Omitted from the current comparison PDF. Stable identifier: `optim-nonsmooth--prox-proj`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id optim-nonsmooth--prox-proj
```

Matching asset directory: `figures/chapter-nonsmooth-convex-optimization/prox-proj`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
