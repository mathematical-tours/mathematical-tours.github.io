# Figure 10.5

Sparse Regularization.

## Figure 10.5

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Sparse spike recovery by the pseudoinverse and by $\lun $ regularization.\relax
```

All six panels arise from one discrete periodic spike experiment. The pseudoinverse discards only numerically zero Fourier modes (tolerance 1e-12), so its noise amplification is visible. Sparse recovery uses the true adjoint and an oracle-selected lambda on a stated grid.

Omitted from the current comparison PDF. Stable identifier: `sparse-regularization--spikes-lun`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id sparse-regularization--spikes-lun
```

Matching asset directory: `figures/chapter-sparse-regularization/spikes-lun`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
