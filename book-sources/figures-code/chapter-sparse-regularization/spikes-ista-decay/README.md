# Figure 10.6

Sparse Regularization.

## Figure 10.6

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Convergence of the energy and iterates under iterative soft thresholding.\relax
```

Actual ISTA objective values and iterates are compared with a much longer accelerated solve of the same objective. The reference is labeled numerical rather than exact. The unit step is valid because max|H|²=1; monotonicity and the reference proximal fixed-point residual are recorded.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--spikes-ista-decay`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-regularization--spikes-ista-decay
```

Matching asset directory: `figures/chapter-sparse-regularization/spikes-ista-decay`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
