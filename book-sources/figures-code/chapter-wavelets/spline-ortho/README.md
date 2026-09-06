# Figure 4.3

Wavelets.

## Figure 4.3

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Orthogonalization of B-splines to define cardinal orthogonal spline functions.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The formula identifies the handwritten example as the orthogonal linear spline. Its triangular generator has compact support, while the orthogonalized function has alternating exponentially decaying tails. The replacement uses the correct Fourier multiplier sqrt(3/(2+cos omega)); its cosine coefficients were evaluated by 8192-point midpoint quadrature. Axis ticks distinguish the generator interval [-1,1] from the displayed orthogonalized interval [-4,4], which uses a compressed horizontal scale. Four small degree examples recover the original spline-family sketches.

Omitted from the current comparison PDF. Stable identifier: `wavelets--spline-ortho`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--spline-ortho
```

Matching asset directory: `figures/chapter-wavelets/spline-ortho`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
