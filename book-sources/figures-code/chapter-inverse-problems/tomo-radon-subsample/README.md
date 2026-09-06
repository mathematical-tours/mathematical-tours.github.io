# Figure 9.8

Inverse Problems.

## Figure 9.8

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Tomography test phantom, discrete Radon projections, and a Cartesian radial Fourier-sampling mask.\relax
```

The middle panel is an actual discrete Radon transform. The last panel illustrates the separate Cartesian radial Fourier model used for reconstruction in Figures 9.9–9.10. The Fourier slice theorem motivates this model, but rasterized radial sampling is not asserted to equal the interpolation-based discrete Radon transform exactly.

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--tomo-radon-subsample`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id inverse-problems--tomo-radon-subsample
```

Matching asset directory: `figures/chapter-inverse-problems/tomo-radon-subsample`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
