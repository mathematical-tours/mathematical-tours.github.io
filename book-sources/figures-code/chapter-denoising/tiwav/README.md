# Figure 7.18

Denoising.

## Figure 7.18

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.19**.

Exact current book caption (LaTeX):

```tex
Translation invariant wavelet coefficients.\relax
```

Shows the actual undecimated two-dimensional Daubechies-2 transform. Every detail array has the original image shape; no subsampling is performed. H means high-pass in the first array coordinate (key da), V in the second (ad), and D in both (dd), matching the tensor-factor convention. Detail grays are signed; the coarse array is divided by four for display.

Omitted from the current comparison PDF. Stable identifier: `denoising--tiwav`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--tiwav
```

Matching asset directory: `figures/chapter-denoising/tiwav`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
