# Figure 4.27

Wavelets.

## Figure 4.27

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.26**.

Exact current book caption (LaTeX):

```tex
Location of large coefficients for Haar (left) and Daubechies-4 (right) wavelets applied to the same piecewise-smooth signal.\relax
```

Haar and Daubechies-4 coefficients of the same piecewise-smooth signal are compared at two common detail levels. Dashed lines identify the input jumps; the coefficient grid uses the same periodic packing convention. Daubechies support spans more samples and is not identified with the one-cell Haar support. Rows share vertical ranges.

Omitted from the current comparison PDF. Stable identifier: `wavelets--vanishing-moments-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--vanishing-moments-1d
```

Matching asset directory: `figures/chapter-wavelets/vanishing-moments-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
