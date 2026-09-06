# Figure 4.19

Wavelets.

## Figure 4.19

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.18**.

Exact current book caption (LaTeX):

```tex
2-D Haar approximation $P_{V_j^O}f$ for increasing $j$.\relax
```

Every displayed image is a true Haar orthogonal projection obtained by averaging the supplied flower image on dyadic squares. Increasing j enlarges the squares and reduces the number of degrees of freedom; all four panels share the same intensity range.

Omitted from the current comparison PDF. Stable identifier: `wavelets--multires-2d-haar`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--multires-2d-haar
```

Matching asset directory: `figures/chapter-wavelets/multires-2d-haar`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
