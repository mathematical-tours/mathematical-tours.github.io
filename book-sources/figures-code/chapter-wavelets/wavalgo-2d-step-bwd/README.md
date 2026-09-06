# Figure 4.24

Wavelets.

## Figure 4.24

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.23**.

Exact current book caption (LaTeX):

```tex
Backward 2-D filterbank step.\relax
```

The inverse graph reverses the analysis order: first reconstruct along the second coordinate from (a_j,d_j^V) and (d_j^H,d_j^D), then along the first coordinate. Each branch inserts zeros before convolution with the unreversed synthesis filter. Plus signs mark the two intermediate sums and the final sum.

Omitted from the current comparison PDF. Stable identifier: `wavelets--wavalgo-2d-step-bwd`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--wavalgo-2d-step-bwd
```

Matching asset directory: `figures/chapter-wavelets/wavalgo-2d-step-bwd`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
