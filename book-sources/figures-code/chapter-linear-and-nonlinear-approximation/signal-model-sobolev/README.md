# Figure 5.5

Linear and Nonlinear Approximation.

## Figure 5.5

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Progressive Gaussian filtering of one white-noise image. Smoothing widths increase from left to right; the displayed images have the same mean and variance.\relax
```

One seeded pure white-noise image is progressively convolved with wider periodic Gaussian kernels. For comparison, each result is affinely normalized to the same mean and standard deviation. The measured gradient energy decreases with smoothing even after this contrast normalization.

Omitted from the current comparison PDF. Stable identifier: `approximation--signal-model-sobolev`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--signal-model-sobolev
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/signal-model-sobolev`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
