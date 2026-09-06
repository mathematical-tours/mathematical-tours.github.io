# Figure 2.1

Fourier and Convolution.

## Figure 2.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Fourier modes (left) and several translated Daubechies wavelets at each scale (right).\relax
```

Each right panel overlays several integer translations of the same actual Daubechies-4 wavelet at a fixed scale. The displayed atoms are 2^(j/2) psi(2^j t-n); supports stay inside the common unit interval. Distinct colors identify translations, while successive rows change frequency or scale.

Omitted from the current comparison PDF. Stable identifier: `fourier--fourier-wav`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id fourier--fourier-wav
```

Matching asset directory: `figures/chapter-fourier-and-convolution/fourier-wav`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
