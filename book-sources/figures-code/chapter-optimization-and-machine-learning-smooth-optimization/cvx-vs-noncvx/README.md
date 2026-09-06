# Figure 14.4

Optimization \& Machine Learning: Smooth Optimization.

## Figure 14.4

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Convex and nonconvex functions; strictly convex and convex but not strictly convex functions. Here $x_t=(1-t)x+ty$ and $\ell _t=(1-t)f(x)+tf(y)$, with $0<t<1$.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Compact exact secant example for nonconvex. The function is .1*x^4-.7*x^2+2, endpoints (-1.8, 1.8), and intermediate point 0. The chord value is 0.78176 and function value 2.0. The chapter caption defines x_t=(1-t)x+ty and ell_t=(1-t)f(x)+t f(y). The red segment is the secant, and its vertical gap shows the defining inequality or its failure. Compact exact secant example for convex. The function is .2*x^2+.6, endpoints (-1.8, 1.8), and intermediate point 0. The chord value is 1.248 and function value 0.6. Compact exact secant example for strictly convex. The function is .2*x^2+.6, endpoints (-1.8, 1.8), and intermediate point 0.6. The chord value is 1.248 and function value 0.672. Compact exact secant example for convex, not strict. The function is .2*max(abs(x)-1,0)^2+.6, endpoints (-0.9, 0.9), and intermediate point 0.2. The chord value is 0.6 and function value 0.6.

Omitted from the current comparison PDF. Stable identifier: `optim-ml-smooth--cvx-vs-noncvx`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id optim-ml-smooth--cvx-vs-noncvx
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-smooth-optimization/cvx-vs-noncvx`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
