# Figure 7.28

Denoising.

## Figure 7.28

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.29**.

Exact current book caption (LaTeX):

```tex
Exact Poisson variances of the Anscombe, Freeman--Tukey, and $2\sqrt {x}$ transforms as functions of the mean intensity.\relax
```

Computes each variance directly from the Poisson probability mass function, including the transformed expectation rather than replacing it by the transform of the mean. The summation tail is below 1e-14. The horizontal line is the asymptotic unit variance; the plot also resolves the low-count regime where stabilization is imperfect.

Omitted from the current comparison PDF. Stable identifier: `denoising--variance-stabilization`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--variance-stabilization
```

Matching asset directory: `figures/chapter-denoising/variance-stabilization`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
