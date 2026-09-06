# Figure 10.3; Figure 10.3 (continued)

Sparse Regularization.

## Figure 10.3

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Panel 1: the objective $\norm {\cdot -y}^2+T^2 \norm {\cdot }_0$. Panels 2--5: the scalar objective $F(x) \eqdef \frac {1}{2}|x-y|^2+\la |x|$ for increasing $\la $.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Exact hard-penalized objective F(x)=(x-y)^2+T^2 1_(x!=0) with y=1.8,T=1.2. The isolated value F(0)=y^2 lies below the open limiting value y^2+T^2, and the unique minimizer is y since y>T. No vertical segment bridges the discontinuity. In general the minimizer is 0 for |y|<T and y for |y|>T, with both minimizing at equality. Exact scalar objective F(x)=0.5(x-y)^2+lambda|x| with y=2 and lambda=0. The minimizer is max(y-lambda,0). Dashed red segments, when present, are the left and right tangents at zero with slopes -y-lambda and lambda-y. Compact labels retain the regime and minimizer; the common objective is given in the main caption.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--varspars`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-regularization--varspars
```

## Figure 10.3 (continued)

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Panel 1: the objective $\norm {\cdot -y}^2+T^2 \norm {\cdot }_0$. Panels 2--5: the scalar objective $F(x) \eqdef \frac {1}{2}|x-y|^2+\la |x|$ for increasing $\la $.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Exact scalar objective F(x)=0.5(x-y)^2+lambda|x| with y=2 and lambda=0.7. The minimizer is max(y-lambda,0). Dashed red segments, when present, are the left and right tangents at zero with slopes -y-lambda and lambda-y. Compact labels retain the regime and minimizer; the common objective is given in the main caption. Exact scalar objective F(x)=0.5(x-y)^2+lambda|x| with y=2 and lambda=2. Exact scalar objective F(x)=0.5(x-y)^2+lambda|x| with y=2 and lambda=2.8.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--varspars-continued`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-regularization--varspars-continued
```

Matching asset directory: `figures/chapter-sparse-regularization/varspars`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
