# Figure 13.20

Basics of Machine Learning.

## Figure 13.20

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Regression using a Gaussian kernel.\relax
```

Restores the original five-panel two-dimensional example: observed colored samples, then four Gaussian-kernel fits with the original bandwidths 0.1, 0.5, 1 and 5. All panels use the same coordinates and response color scale. Coefficients solve (K+0.03 I)c=y for one seeded noisy sample, and predictions evaluate the same kernel expansion on the grid.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--kernel`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--kernel
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/kernel`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
