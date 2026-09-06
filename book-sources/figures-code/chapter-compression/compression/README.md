# Figure 6.2

Compression.

## Figure 6.2

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Image compression using wavelet support coding.\relax
```

Recomputes the orthogonal wavelet coefficients, dead-zone integer quantization, support, and midpoint reconstruction. The binary example stores a header, unsigned 16-bit indices, and signed 16-bit nonzero symbols; its measured length includes the header. It is an explicit support codec, not JPEG-2000. Zero is reconstructed as zero, and nonzero symbols as sign(q)(|q|+1/2)T.

Omitted from the current comparison PDF. Stable identifier: `compression--compression`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id compression--compression
```

Matching asset directory: `figures/chapter-compression/compression`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
