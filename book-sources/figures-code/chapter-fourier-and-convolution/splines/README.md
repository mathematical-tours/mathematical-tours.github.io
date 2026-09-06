# Figure 2.7

Fourier and Convolution.

## Figure 2.7

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Cardinal splines are defined by successive convolutions.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The box, triangle, and quadratic spline are exact centered cardinal B-splines B0, B1=B0*B0, and B2=B1*B0. Degree and convolution labels specify each spline without ambiguity.

Omitted from the current comparison PDF. Stable identifier: `fourier--splines`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id fourier--splines
```

Matching asset directory: `figures/chapter-fourier-and-convolution/splines`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
