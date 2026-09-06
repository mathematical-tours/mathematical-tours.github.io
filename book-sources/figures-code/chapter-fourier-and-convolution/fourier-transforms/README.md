# Figure 2.2

Fourier and Convolution.

## Figure 2.2

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Four settings for Fourier analysis, linked by sampling and periodization.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Each card pairs a signal domain with its Fourier domain. Spatial periodization acts horizontally and sampling vertically. With h=2pi/N, both paths give b_n=p(hn)=sum_j a_(n+jN). The formulas use the chapter convention: c_k=f-hat(k)/(2pi), A(theta)=(1/h) sum_j f-hat((theta+2pi j)/h), and B_k=A(2pi k/N)=N sum_j c_(k+jN). Rapid decay and smoothness ensure convergence. The unit-period and 2pi-period conventions are consistent. The forward transforms are not unitary with the displayed normalizations.

Omitted from the current comparison PDF. Stable identifier: `fourier--fourier-transforms`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id fourier--fourier-transforms
```

Matching asset directory: `figures/chapter-fourier-and-convolution/fourier-transforms`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
