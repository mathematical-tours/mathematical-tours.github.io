# Figure 4.22

Wavelets.

## Figure 4.22

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.21**.

Exact current book caption (LaTeX):

```tex
Forward 2-D filterbank step.\relax
```

The 2-D analysis graph follows the chapter’s definitions exactly: H is filtering/downsampling in the first coordinate and V in the second. Low-low yields a_j, low-high yields d_j^V, high-low yields d_j^H, and high-high yields d_j^D. Every downsampling follows convolution with the reversed analysis filter.

Omitted from the current comparison PDF. Stable identifier: `wavelets--wavalgo-2d-step-fwd`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--wavalgo-2d-step-fwd
```

Matching asset directory: `figures/chapter-wavelets/wavalgo-2d-step-fwd`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
