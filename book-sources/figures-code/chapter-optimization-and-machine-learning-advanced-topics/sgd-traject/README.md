# Figure 15.5

Optimization \& Machine Learning: Advanced Topics.

## Figure 15.5

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Repeated SGD runs for $k\mapsto x_k\in \RR $. Left: the distribution of iterates at each iteration. Right: individual trajectories.\relax
```

Simulates SGD for E[(x-Z)²/2] with Z standard Gaussian, whose minimizer is exactly zero. Every step draws fresh independent noise. The left heatmap is the empirical density of the same trials whose first 18 trajectories appear at right; the deterministic starting point is represented by one histogram bin; a logarithmic color scale resolves the later spread without allowing this initial point mass to dominate the range.

Omitted from the current comparison PDF. Stable identifier: `optim-ml-advanced--sgd-traject`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id optim-ml-advanced--sgd-traject
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-advanced-topics/sgd-traject`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
