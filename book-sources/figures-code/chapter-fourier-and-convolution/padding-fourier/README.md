# Figure 2.12

Fourier and Convolution.

## Figure 2.12

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Interpolation by zero-padding in the frequency domain.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The retained signed-frequency coefficients are padded with zeros, and the inverse transform is multiplied by Q/N. The illustrative trigonometric polynomial and coefficient magnitudes are consistent for N=5 and Q=15; dark markers identify the five original samples and smaller colored markers the finer grid. These are analytic example values, not measurements. Odd N avoids the even-length Nyquist splitting convention.

Omitted from the current comparison PDF. Stable identifier: `fourier--padding-fourier`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id fourier--padding-fourier
```

Matching asset directory: `figures/chapter-fourier-and-convolution/padding-fourier`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
