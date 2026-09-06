# Figure 4.15

Wavelets.

## Figure 4.15

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.14**.

Exact current book caption (LaTeX):

```tex
Inverse filter-bank reconstruction.\relax
```

The synthesis stage inserts zeros before convolution with h and g and sums the branches. These filters are not reversed in synthesis. The formula is the exact adjoint of the forward filtering/downsampling step for the real orthogonal filter convention used in the chapter.

Omitted from the current comparison PDF. Stable identifier: `wavelets--filterbank-bwd`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--filterbank-bwd
```

Matching asset directory: `figures/chapter-wavelets/filterbank-bwd`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
