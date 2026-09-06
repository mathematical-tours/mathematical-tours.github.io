# Figure 4.9

Wavelets.

## Figure 4.9

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Haar filter weights as inner products.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The original inner products are +1 for the scaling function and +1 or -1 for the wavelet. The chapter filter definition adds a factor 1/sqrt(2); the resulting h=(1,1)/sqrt(2) and g=(1,-1)/sqrt(2) are therefore explicitly distinguished from the raw overlap integrals.

Omitted from the current comparison PDF. Stable identifier: `wavelets--haar-inner`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--haar-inner
```

Matching asset directory: `figures/chapter-wavelets/haar-inner`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
