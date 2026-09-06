# Figure 5.9

Linear and Nonlinear Approximation.

## Figure 5.9

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Wavelet approximation error for the images in Figure \ref {fig-approx-speed-1}.\relax
```

Curves are the exact discarded coefficient energies in the same orthonormal Daubechies-4 transform of the four images in Figure 5.8. The vertical view stops at relative squared error 1e-6 to emphasize the approximation regime instead of roundoff; all raw energies remain in numerical-data.npz. This finite experiment is not a fit to asymptotic theoretical exponents.

Omitted from the current comparison PDF. Stable identifier: `approximation--approx-speed`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--approx-speed
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/approx-speed`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
