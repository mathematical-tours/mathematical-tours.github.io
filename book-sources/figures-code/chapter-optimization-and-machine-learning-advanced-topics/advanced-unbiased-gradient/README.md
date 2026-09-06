# Figure 15.4

Optimization \& Machine Learning: Advanced Topics.

## Figure 15.4

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Unbiased stochastic gradient estimate (panel 1) and a schematic SGD trajectory (panel 2).\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The red dashed vector is the exact arithmetic mean of the three individual gradient vectors. Its endpoint is the centroid of the lightly shaded triangle joining their endpoints. Uniform sampling from {1,2,3} is stated explicitly, which is required for this unbiased estimate of the arithmetic-mean objective. The path and shrinking ellipses preserve the qualitative intent of the sketch. They are illustrative neighborhoods, not confidence regions or a simulated dataset. A dotted continuation separates the last displayed iterate from the target minimizer. Convergence requires the objective, noise, and step-size assumptions developed in the section; step sizes alone do not guarantee it.

Omitted from the current comparison PDF. Stable identifier: `optim-ml-advanced--advanced-unbiased-gradient`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id optim-ml-advanced--advanced-unbiased-gradient
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-advanced-topics/advanced-unbiased-gradient`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
