# Figure 9.2

Inverse Problems.

## Figure 9.2

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Left: capping reciprocal singular values. Right: the Tikhonov bound $\mu _{\la }(\si ) \leq C_\la = \frac {1}{2\sqrt {\la }}$.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The left plateau is interpreted as the capped inverse 1/max(sigma,tau), rather than truncated SVD. A separate cutoff tau avoids confusing its parameter with the squared-scale Tikhonov parameter lambda. Both plots now use common axis scales and tau=sqrt(lambda): the Tikhonov maximum at sigma=sqrt(lambda) is half the capped gain. The dashed reciprocal curves have the same scale.

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--bound-regul-variance`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id inverse-problems--bound-regul-variance
```

Matching asset directory: `figures/chapter-inverse-problems/bound-regul-variance`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
