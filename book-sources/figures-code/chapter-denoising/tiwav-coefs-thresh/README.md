# Figure 7.19

Denoising.

## Figure 7.19

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.20**.

Exact current book caption (LaTeX):

```tex
Left: translation invariant wavelet coefficients, for $j=-8, \om =H$, right: thresholded coefficients.\relax
```

Uses a 512 by 512 flower image so the first detail level is j=-8 in the chapter’s J=-log2(N) convention. The H band is high-pass in the first array coordinate. The same signed color range [-4 sigma,4 sigma] is used, with saturation only for display, before and after exact hard thresholding; equality at the threshold maps to zero.

Omitted from the current comparison PDF. Stable identifier: `denoising--tiwav-coefs-thresh`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--tiwav-coefs-thresh
```

Matching asset directory: `figures/chapter-denoising/tiwav-coefs-thresh`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
