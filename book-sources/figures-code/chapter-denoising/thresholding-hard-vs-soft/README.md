# Figure 7.11

Denoising.

## Figure 7.11

Previously accepted TikZ drawing retained in the reading editions.

Previous audit identifier: **7.12**.

Exact current book caption (LaTeX):

```tex
Hard and soft thresholding functions. \relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Solid hard thresholding has open diagonal endpoints and closed zero endpoints at +/-T. Dashed soft thresholding is continuous and subtracts T from the retained magnitude. No vertical line is drawn across a discontinuity; both functions are zero exactly at the cutoffs.

Omitted from the current comparison PDF. Stable identifier: `denoising--thresholding-hard-vs-soft`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--thresholding-hard-vs-soft
```

Matching asset directory: `figures/chapter-denoising/thresholding-hard-vs-soft`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
