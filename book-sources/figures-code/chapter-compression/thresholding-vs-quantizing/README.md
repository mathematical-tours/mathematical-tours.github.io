# Figure 6.1

Compression.

## Figure 6.1

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Hard thresholding and quantization as functions of the input coefficient.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The horizontal coordinate is x/T; the hard-threshold output is divided by T so both vertical axes are dimensionless. Hard thresholding sends +/-T to zero, whereas the chapter integer quantizer sends them to +/-1. All open and closed endpoints follow these definitions. In particular the quantizer zero bin is open at both ends.

Omitted from the current comparison PDF. Stable identifier: `compression--thresholding-vs-quantizing`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id compression--thresholding-vs-quantizing
```

Matching asset directory: `figures/chapter-compression/thresholding-vs-quantizing`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
