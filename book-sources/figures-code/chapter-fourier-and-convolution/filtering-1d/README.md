# Figure 2.5

Fourier and Convolution.

## Figure 2.5

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Filtering an irregular signal with Gaussian filters of increasing width $\si $.\relax
```

The three smoothed curves are computed from one piecewise smooth signal by convolution with normalized discrete Gaussian kernels, using periodic boundary conditions. The displayed standard deviations are in physical time units. The original signal is retained as a thin guide so widening the filter visibly trades detail for smoothness.

Omitted from the current comparison PDF. Stable identifier: `fourier--filtering-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id fourier--filtering-1d
```

Matching asset directory: `figures/chapter-fourier-and-convolution/filtering-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
