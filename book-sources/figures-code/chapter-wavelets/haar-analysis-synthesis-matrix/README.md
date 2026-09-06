# Figure 4.14

Wavelets.

## Figure 4.14

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.13 (additional)**.

Exact current book caption (LaTeX):

```tex
One-step Haar analysis matrix $W$ and synthesis matrix $W^{\mathsf T}$ for eight samples.\relax
```

Displays the actual N=8 one-step Haar analysis matrix and its transpose. Nonzero entries are signed 1/sqrt(2), and W^T W=I is checked. The layout makes explicit that an orthogonal inverse is the transpose and reverses the order of composed stages.

Omitted from the current comparison PDF. Stable identifier: `wavelets--haar-analysis-synthesis-matrix`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--haar-analysis-synthesis-matrix
```

Matching asset directory: `figures/chapter-wavelets/haar-analysis-synthesis-matrix`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
