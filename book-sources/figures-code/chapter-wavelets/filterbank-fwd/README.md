# Figure 4.10

Wavelets.

## Figure 4.10

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Forward filter bank decomposition.\relax
```

The two analysis branches convolve with reversed real filters h-bar and g-bar before retaining even samples. Only the approximation branch is decomposed again. Indices match the chapter: increasing j is coarser. The displayed coefficient formula fixes the convolution and downsampling convention.

Omitted from the current comparison PDF. Stable identifier: `wavelets--filterbank-fwd`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--filterbank-fwd
```

Matching asset directory: `figures/chapter-wavelets/filterbank-fwd`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
