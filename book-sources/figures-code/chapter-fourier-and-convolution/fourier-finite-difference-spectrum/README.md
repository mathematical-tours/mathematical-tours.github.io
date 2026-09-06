# Figure 2.19

Fourier and Convolution.

## Figure 2.19

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Comparison of the spectra of $\De $ and $D_2$.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The positive curves are explicitly labeled as negative eigenvalues divided by N squared. The continuum curve is (2pi k/N) squared, reaching pi squared at Nyquist; the discrete curve is 4 sin squared(pi k/N), reaching four. Actual Laplacian eigenvalues are nonpositive.

Omitted from the current comparison PDF. Stable identifier: `fourier--fourier-finite-difference-spectrum`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id fourier--fourier-finite-difference-spectrum
```

Matching asset directory: `figures/chapter-fourier-and-convolution/fourier-finite-difference-spectrum`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
