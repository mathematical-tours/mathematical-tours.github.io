# Figure 13.1

Basics of Machine Learning.

## Figure 13.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Empirical covariance of the Iris measurements and its eigenvalues.\relax
```

The 150 real Iris observations are centered without feature standardization. Covariance uses the chapter normalization 1/n. The spectrum displays covariance eigenvalues sigma_k²/n for singular values of the unnormalized centered data, avoiding a missing factor n.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--cov`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--cov
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/cov`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
