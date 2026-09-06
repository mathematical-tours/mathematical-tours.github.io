# Figure 5.23

Linear and Nonlinear Approximation.

## Figure 5.23

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
A curvelet $c_m$ in space (left) and its localization in the Fourier domain (right).\relax
```

The real atom is generated from a smooth compact annular wedge with parabolic radial/transverse bandwidths, and its Fourier image is computed from that same atom. Taking the real part introduces the conjugate wedge. It belongs to the explicit periodic oversampled curvelet construction used for Figure 5.26; it is not an arbitrary Gaussian labeled a tight-frame curvelet.

Omitted from the current comparison PDF. Stable identifier: `approximation--curvelet`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--curvelet
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/curvelet`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
