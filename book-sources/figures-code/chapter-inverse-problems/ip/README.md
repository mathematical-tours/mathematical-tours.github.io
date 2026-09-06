# Figure 9.1

Inverse Problems.

## Figure 9.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Examples of forward operators in inverse problems.\relax
```

Forward operators are applied to a common acquired image: periodic Gaussian convolution, restriction to one pixel in each 4 by 4 cell, and a known random observation mask. These are measurements, not reconstructions.

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--ip`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id inverse-problems--ip
```

Matching asset directory: `figures/chapter-inverse-problems/ip`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
