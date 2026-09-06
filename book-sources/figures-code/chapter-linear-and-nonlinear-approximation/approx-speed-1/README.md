# Figure 5.8

Linear and Nonlinear Approximation.

## Figure 5.8

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Test images with different structures: smooth, cartoon, flower, and mandrill.\relax
```

These four precisely defined test images are reused unchanged in the next approximation-error figure: smooth, cartoon, the supplied flower and the author-requested mandrill.

Omitted from the current comparison PDF. Stable identifier: `approximation--approx-speed-1`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--approx-speed-1
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/approx-speed-1`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
