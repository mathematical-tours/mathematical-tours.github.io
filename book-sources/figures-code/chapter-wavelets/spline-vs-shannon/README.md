# Figure 4.6

Wavelets.

## Figure 4.6

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Spline and Shannon wavelet segmentation of the frequency axis.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The Shannon row shows the exact low-pass band |omega|<=pi and the next two dyadic wavelet bands. The spline row deliberately shows schematic smooth overlap, not the spectrum of a named spline order. Only a finite frequency window is drawn; the partition identity refers to all scales. Cropped neighboring prose in the scan is omitted. The energy of a scale-j wavelet is multiplied by 2 to the power -j to remove its L2 dilation factor, as the displayed identity explains; this gives the factor two on the outer scale-minus-one band.

Omitted from the current comparison PDF. Stable identifier: `wavelets--spline-vs-shannon`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--spline-vs-shannon
```

Matching asset directory: `figures/chapter-wavelets/spline-vs-shannon`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
