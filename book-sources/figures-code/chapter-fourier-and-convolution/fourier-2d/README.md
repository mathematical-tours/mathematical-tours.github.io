# Figure 2.17

Fourier and Convolution.

## Figure 2.17

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
2-D Fourier analysis of an image (left), and attenuation of the periodicity artifact using masking (right).\relax
```

Both centered spectra are recomputed from the displayed flower images and share a decibel reference and color scale. A separable Hann mask forces the image toward zero on the boundary, reducing discontinuities in its periodic extension. Windowing also changes the image and broadens spectral features; it is not a lossless correction.

Omitted from the current comparison PDF. Stable identifier: `fourier--fourier-2d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id fourier--fourier-2d
```

Matching asset directory: `figures/chapter-fourier-and-convolution/fourier-2d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
