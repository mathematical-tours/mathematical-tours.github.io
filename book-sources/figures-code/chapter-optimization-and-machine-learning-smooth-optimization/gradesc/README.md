# Figure 14.10

Optimization \& Machine Learning: Smooth Optimization.

## Figure 14.10

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Effect of the step size $\tau $ on gradient descent (panel 1) and exact line search (panel 2).\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Both trajectories are exact gradient-descent iterates for the displayed quadratic. The larger step oscillates but remains stable because both step sizes lie in (0,2/L), with L=5. The sketch is schematic; no numerical values were legible in it. The path uses exact line search for the displayed positive definite quadratic f(x)=(x₁ squared+5x₂ squared)/2. The right-angle marker is computed from consecutive actual steps. The later, shorter steps remain connected as the path approaches the minimizer.

Omitted from the current comparison PDF. Stable identifier: `optim-ml-smooth--gradesc`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id optim-ml-smooth--gradesc
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-smooth-optimization/gradesc`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
