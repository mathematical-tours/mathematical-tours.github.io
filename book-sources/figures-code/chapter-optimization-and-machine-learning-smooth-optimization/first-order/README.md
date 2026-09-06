# Figure 14.7

Optimization \& Machine Learning: Smooth Optimization.

## Figure 14.7

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Local extrema (panel 1), a stationary inflection point (panel 2), and a global minimizer (panel 3).\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The exact quartic f(x)=(x^2-1)^2+0.65 has local minima at -1 and 1 and a local maximum at 0. Horizontal tangent segments identify the three stationary points. The exact cubic f(x)=x^3 has a horizontal tangent at zero and a stationary inflection, which is neither a minimum nor a maximum. Plot coordinates have a vertical drawing offset, with the horizontal axis at f=0. The tangent plane touches the bowl at its minimum and supports the graph from below. The plane is horizontal in the depicted three-dimensional coordinates.

Omitted from the current comparison PDF. Stable identifier: `optim-ml-smooth--first-order`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id optim-ml-smooth--first-order
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-smooth-optimization/first-order`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
