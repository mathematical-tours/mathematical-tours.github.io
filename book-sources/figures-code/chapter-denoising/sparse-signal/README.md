# Figure 7.14

Denoising.

## Figure 7.14

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.15**.

Exact current book caption (LaTeX):

```tex
Left: clean sparse coefficients $a_0$; right: noisy coefficients $a$.\relax
```

A 256-dimensional coefficient vector with exactly twelve nonzeros is corrupted by independent Gaussian noise of standard deviation 0.12. Signal and noise are generated from a fixed seed; both plots have the same axes.

Omitted from the current comparison PDF. Stable identifier: `denoising--sparse-signal`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--sparse-signal
```

Matching asset directory: `figures/chapter-denoising/sparse-signal`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
