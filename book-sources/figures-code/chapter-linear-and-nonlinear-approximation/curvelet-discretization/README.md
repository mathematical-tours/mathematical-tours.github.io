# Figure 5.24

Linear and Nonlinear Approximation.

## Figure 5.24

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Sampling pattern for the curvelet positions.\relax
```

Centers follow exactly u=R_theta(2^(j/2)m1,2^j m2), translated to the middle of the plotting window. The highlighted support has width equal to length squared. Refining scale changes the two sampling intervals differently; rotation rotates the whole lattice.

Omitted from the current comparison PDF. Stable identifier: `approximation--curvelet-discretization`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--curvelet-discretization
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/curvelet-discretization`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
