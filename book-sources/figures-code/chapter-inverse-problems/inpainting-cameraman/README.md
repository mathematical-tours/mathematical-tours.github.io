# Figure 9.6

Inverse Problems.

## Figure 9.6

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Inpainting with Sobolev and TV regularization.\relax
```

Sobolev and isotropic nonsmooth-TV inpainting share the exact same random observation mask. Known pixels are imposed as equality constraints at every iteration. Recomputed SNR values determine the outcome rather than preserving the ranking of a different legacy dataset.

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--inpainting-cameraman`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id inverse-problems--inpainting-cameraman
```

Matching asset directory: `figures/chapter-inverse-problems/inpainting-cameraman`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
