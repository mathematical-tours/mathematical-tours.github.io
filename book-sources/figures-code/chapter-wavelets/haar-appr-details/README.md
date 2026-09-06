# Figure 4.5

Wavelets.

## Figure 4.5

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Three successive Haar approximations (top) and their extracted details (bottom). Coarser approximations use fewer intervals.\relax
```

Three successive Haar resolution pairs are arranged above their corresponding detail functions, including the additional 64-to-32 interval pair. Each column satisfies finer = coarser + detail and exact coarse/detail orthogonality. The resolution and detail rows use common vertical ranges.

Omitted from the current comparison PDF. Stable identifier: `wavelets--haar-appr-details`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--haar-appr-details
```

Matching asset directory: `figures/chapter-wavelets/haar-appr-details`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
