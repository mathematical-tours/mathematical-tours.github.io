# Figure 2.18

Fourier and Convolution.

## Figure 2.18

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Heat diffusion as a convolution.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The left curves show an illustrative signal with two Fourier modes at time zero and after positive diffusion times. The right curves are Gaussian heat kernels for the same two times: increasing time broadens the kernel and lowers its peak. The normalization and variance follow the chapter heat equation.

Omitted from the current comparison PDF. Stable identifier: `fourier--heat`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id fourier--heat
```

Matching asset directory: `figures/chapter-fourier-and-convolution/heat`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
