# Figure 10.4

Sparse Regularization.

## Figure 10.4

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Geometry of convex optimization: smooth-boundary and corner contact with objective level sets, followed by the smallest scaled $\ell ^2$ and $\ell ^1$ balls touching an affine constraint. Red dots mark the contact points.\relax
```

The first two panels show exact constraint contact points for F(x)=(x1-1.8)^2+2(x2-.3)^2. At the l1 corner, the supporting outward normal lies in the diamond normal cone. For x1+0.6x2=1, the minimal Euclidean ball touches at (1,.6)/1.36 and the minimal l1 ball at (1,0); both contacts and their correctly scaled balls are displayed.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--optim-geom`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id sparse-regularization--optim-geom
```

Matching asset directory: `figures/chapter-sparse-regularization/optim-geom`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
