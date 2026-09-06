# Figure 4.12

Wavelets.

## Figure 4.12

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Packed coefficient arrays after one, two, and three wavelet decomposition steps.\relax
```

The storage arrays are recomputed after one, two and three Haar steps. The active approximation occupies the leftmost block and is replaced by a smaller approximation followed by its detail. Earlier details stay in place. Length and energy are preserved at every step.

Omitted from the current comparison PDF. Stable identifier: `wavelets--wavalgo-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--wavalgo-1d
```

Matching asset directory: `figures/chapter-wavelets/wavalgo-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
