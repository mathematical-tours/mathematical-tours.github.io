# Figure 5.6

Linear and Nonlinear Approximation.

## Figure 5.6

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Progressive dead-leaves construction: later opaque shapes cover earlier ones, producing a family of cartoon images.\relax
```

A progressive dead-leaves model overlays opaque ellipses with independently sampled positions, sizes, orientations and gray values. Later leaves occlude earlier ones. Every snapshot is piecewise constant with a finite union of smooth boundary arcs; intersections and occlusions create corners. Total variation is measured, not assumed monotone under occlusion.

Omitted from the current comparison PDF. Stable identifier: `approximation--image-model-cartoon`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--image-model-cartoon
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/image-model-cartoon`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
