# Figure 4.7

Wavelets.

## Figure 4.7

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Periodized basis functions. Only the translations $0\leq n<2^{-j}$ are distinct, and $j_0=0$ is a common choice of coarsest scale.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The illustrative compact bump is cos squared(2pi(x-0.87)) for |x-0.87|<=0.25, zero elsewhere. On [0,1], only f(x) and its left-shifted copy f(x+1) contribute, so the drawn sum and matching endpoint values are exact. The profile is not presented as a particular orthogonal wavelet; periodization applies to each basis function individually. No unknown profile or empirical values are inferred from the scan. The selected translation centers lie in the half-open unit interval. The example j=-2 has four distinct centers; the center at one is the same periodic translate as zero. Smooth pulse profiles are schematic and do not assert that this particular shape generates an orthogonal wavelet basis.

Omitted from the current comparison PDF. Stable identifier: `wavelets--periodize`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--periodize
```

Matching asset directory: `figures/chapter-wavelets/periodize`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
