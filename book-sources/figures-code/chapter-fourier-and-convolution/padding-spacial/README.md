# Figure 2.11

Fourier and Convolution.

## Figure 2.11

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Fourier transform approximation by zero-padding in the spatial domain.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. An exact analytic example uses N=12, Q=24 and f(x)=sin squared(pi x) on [0,1], zero elsewhere. The Q spatial samples stop at T-1/N, where T=Q/N; the final point T is excluded. The continuous curve is the magnitude of the actual Fourier transform, and the Q markers are normalized DFT magnitudes at signed indices -Q/2 through Q/2-1. Thus the unchanged Nyquist bound is pi N and the frequency step is 2pi/T, correcting the handwritten reciprocal bounds. The DFT approximates the continuous transform by a Riemann sum; it does not equal it exactly. These are analytic values, not measurements.

Omitted from the current comparison PDF. Stable identifier: `fourier--padding-spacial`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id fourier--padding-spacial
```

Matching asset directory: `figures/chapter-fourier-and-convolution/padding-spacial`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
