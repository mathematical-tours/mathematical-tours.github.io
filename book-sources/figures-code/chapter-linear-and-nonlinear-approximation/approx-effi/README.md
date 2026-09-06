# Figure 5.11

Linear and Nonlinear Approximation.

## Figure 5.11

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Comparison of approximation error decay for different bases.\relax
```

Actual best-M errors in four orthonormal bases use one common image and real coefficient counts. The local DCT consists of disjoint 16 by 16 blocks. Full-transform inversion is checked before computing the discarded-energy curves.

Omitted from the current comparison PDF. Stable identifier: `approximation--approx-effi`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--approx-effi
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/approx-effi`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
