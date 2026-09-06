# Figure 5.12

Linear and Nonlinear Approximation.

## Figure 5.12

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Summary of linear and nonlinear approximation rates $\norm {f-f_M}^2$ for different classes of 1-D signals and images.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The sketches are reorganized into an upper-rate table. L and NL mean linear and nonlinear approximation. The Sobolev columns explicitly refer to H^alpha, avoiding a conflation with noninteger Holder regularity at the endpoint. The one-dimensional piecewise-smooth rates are shown for alpha at least one; wavelets require the stated regularity, moments, and boundary assumptions. The BV class also has the chapter uniform bound. The cartoon column assumes C2 pieces and contours, illustrated by a smooth closed curve. The curvelet rate retains the logarithmic factor omitted by the handwritten summary. The entries are upper rates, not assertions that every signal attains them.

Omitted from the current comparison PDF. Stable identifier: `approximation--approx-comparaison`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--approx-comparaison
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/approx-comparaison`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
