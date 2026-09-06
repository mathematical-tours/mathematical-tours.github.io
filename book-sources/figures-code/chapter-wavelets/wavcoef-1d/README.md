# Figure 4.8

Wavelets.

## Figure 4.8

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Wavelet coefficients. The first column shows the input signal and its packed decomposition; the other panels show individual detail scales.\relax
```

The first column contains the input and complete packed Haar transform. The next two columns enlarge its four detail bands, with matching scale colors. All coefficients come from the same orthonormal transform; coarse coefficients precede details from coarse to fine.

Omitted from the current comparison PDF. Stable identifier: `wavelets--wavcoef-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--wavcoef-1d
```

Matching asset directory: `figures/chapter-wavelets/wavcoef-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
