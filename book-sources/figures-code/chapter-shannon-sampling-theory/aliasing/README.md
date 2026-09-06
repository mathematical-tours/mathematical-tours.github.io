# Figure 1.7

Shannon Sampling Theory.

## Figure 1.7

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Aliasing of a sine wave. This periodic example is interpreted through its discrete spectral lines rather than as an $L^2(\RR )$ signal.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. For unit sample spacing, cos(pi x/2) and cos(3 pi x/2) agree at every integer. Frequencies are identified modulo 2 pi; the high-frequency pair folds to the low-frequency pair. This example is periodic and is not an L2 signal on the real line.

Omitted from the current comparison PDF. Stable identifier: `shannon--aliasing`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id shannon--aliasing
```

Matching asset directory: `figures/chapter-shannon-sampling-theory/aliasing`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
