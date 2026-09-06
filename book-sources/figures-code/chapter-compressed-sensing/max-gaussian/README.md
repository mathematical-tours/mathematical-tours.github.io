# Figure 12.4

Compressed Sensing.

## Figure 12.4

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Maximum of $N$ independent standard Gaussian variables, compared with the dashed curve $\sqrt {2\log (N)}$.\relax
```

Uses the signed maximum max(Z1,...,ZN), as the caption states, rather than the maximum absolute value shown earlier in the denoising chapter. Independent trials determine the empirical mean and Monte Carlo standard error; the dashed expression is an asymptotic comparison, not an exact finite-N mean.

Omitted from the current comparison PDF. Stable identifier: `compressed-sensing--max-gaussian`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id compressed-sensing--max-gaussian
```

Matching asset directory: `figures/chapter-compressed-sensing/max-gaussian`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
