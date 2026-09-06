# Figure 12.5

Compressed Sensing.

## Figure 12.5

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Phase transitions. Right: empirical probabilities that recovery criteria hold as sparsity varies. Blue: weak ERC; black: ERC; green: $\norm {\eta _F}_\infty \leq 1$; red: exact recovery.\relax
```

All probabilities are Monte Carlo counts from column-normalized Gaussian matrices and independent random supports/signs. Exact recovery is solved by linear programming with relative-scale tolerance 1e-5. ERC tests max_j||A_I^+ a_j||1<1. Weak ERC uses the Neumann-series bound max_j||A_Iᵀa_j||1 / (1-||A_IᵀA_I-I||1)<1, requiring a positive denominator. Fuchs tests the sign-specific off-support certificate. Finite trial counts are stated; these are not asymptotic threshold curves.

Omitted from the current comparison PDF. Stable identifier: `compressed-sensing--phase-trans`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id compressed-sensing--phase-trans
```

Matching asset directory: `figures/chapter-compressed-sensing/phase-trans`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
