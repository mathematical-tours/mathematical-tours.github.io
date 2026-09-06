# Figure 2.4

Fourier and Convolution.

## Figure 2.4

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Signal filtering with a box filter (running average). \relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The irregular illustrative signal is a finite trigonometric sum. Each Fourier mode is integrated exactly over the unit-width box, so the dashed curve is its actual convolution and the marked value is the corresponding integral. The box has mass one.

Omitted from the current comparison PDF. Stable identifier: `fourier--filter-box`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id fourier--filter-box
```

Matching asset directory: `figures/chapter-fourier-and-convolution/filter-box`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
