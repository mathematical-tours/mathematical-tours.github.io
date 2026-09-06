# Figure 5.25

Linear and Nonlinear Approximation.

## Figure 5.25

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Wavelet (left) and curvelet (right) approximation of a cartoon image.\relax
```

Reconstructs the original geometric schematic rather than replacing it with an unrelated numerical panel. Wavelet footprints have comparable tangential/normal scales and a few fixed orientations; curvelet footprints are elongated and follow the contour tangent. Red transverse atoms illustrate an orientation mismatch. Footprints indicate localization, not exact compact supports or numerical coefficient values.

Omitted from the current comparison PDF. Stable identifier: `approximation--curvelets-vs-wavelets`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--curvelets-vs-wavelets
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/curvelets-vs-wavelets`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
