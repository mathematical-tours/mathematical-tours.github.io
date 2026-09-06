# Figure 4.26

Wavelets.

## Figure 4.26

Previously accepted TikZ drawing retained in the reading editions.

Previous audit identifier: **4.25**.

Exact current book caption (LaTeX):

```tex
Vanishing moments correspond to vanishing derivatives of the Fourier transform.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The drawing separates the low-pass zero at pi from the high-pass zero at zero, avoiding the ambiguous overlaid labels in the scan. The curves are magnitudes for a two-moment QMF example; the derivative statement concerns the smooth filter symbols themselves, including their phases. The displayed zeros have order two, while the formula states the general p-moment equivalence.

Omitted from the current comparison PDF. Stable identifier: `wavelets--vanmoments`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--vanmoments
```

Matching asset directory: `figures/chapter-wavelets/vanmoments`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
