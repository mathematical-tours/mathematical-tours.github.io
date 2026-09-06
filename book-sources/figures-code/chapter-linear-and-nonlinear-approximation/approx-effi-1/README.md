# Figure 5.10

Linear and Nonlinear Approximation.

## Figure 5.10

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Fourier, DCT, and wavelet approximations of the same $256\times 256$ flower image, each retaining $M=4096$ real orthonormal coefficients.\relax
```

All three orthonormal representations retain exactly 4096 real coefficients from the same 256 by 256 flower image. Fourier counts individual real sine/cosine atoms; wavelets use four periodic db4 levels. At this larger common budget, wavelets have the lowest measured error. This observation concerns this image and budget, not a universal ordering of bases.

Omitted from the current comparison PDF. Stable identifier: `approximation--approx-effi-1`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--approx-effi-1
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/approx-effi-1`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
