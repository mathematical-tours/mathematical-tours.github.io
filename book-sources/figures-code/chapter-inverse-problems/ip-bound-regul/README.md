# Figure 9.3

Inverse Problems.

## Figure 9.3

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Bounding $\la \frac {\si ^\be }{\la +\si ^2} \leq D_{\la ,\be }$.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Corrects the handwritten maximizer: sigma*=sqrt(beta*lambda/(2-beta)), including the missing factor beta. For beta=2 the supremum lambda is approached asymptotically. For beta>2 the supremum over all sigma>=0 is infinite. The figure explicitly distinguishes this from the bounded spectrum of a fixed bounded operator. Curves show representative beta values with rescaled axes and lambda>0.

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--ip-bound-regul`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id inverse-problems--ip-bound-regul
```

Matching asset directory: `figures/chapter-inverse-problems/ip-bound-regul`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
