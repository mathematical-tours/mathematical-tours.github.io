# Figure 8.5

Variational Priors and Regularization.

## Figure 8.5

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Regularized absolute value $x \mapsto \sqrt {x^2+\epsilon ^2}$.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Exact graphs of |x| and sqrt(x^2+epsilon^2), with epsilon=0.1. Both panels use identical axes near the origin, making the smoothing scale comparable. The smoothed function lies above absolute value and has value epsilon and derivative zero at the origin. The old plots were too coarse to reliably convey the stated epsilon values. Exact graphs of |x| and sqrt(x^2+epsilon^2), with epsilon=0.01.

Omitted from the current comparison PDF. Stable identifier: `variational-priors--regul-abs`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id variational-priors--regul-abs
```

Matching asset directory: `figures/chapter-variational-priors-and-regularization/regul-abs`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
