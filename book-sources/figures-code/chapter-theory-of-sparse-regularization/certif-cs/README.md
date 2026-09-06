# Figure 11.9

Theory of Sparse Regularization.

## Figure 11.9

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Fuchs precertificate $\eta _F$ for a Gaussian matrix $A \in \RR ^{P \times N}$ with $N=64$ and increasing measurement count $P$.\relax
```

The Fuchs precertificate is recomputed by the exact Gram solve Aᵀ A_I(A_Iᵀ A_I)^(-1)sign(x_I). The same four signed support entries and a nested Gaussian array are used as P increases. Dashed ±1 lines show dual feasibility; no off-support clipping is applied.

Omitted from the current comparison PDF. Stable identifier: `sparse-theory--certif-cs`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-theory--certif-cs
```

Matching asset directory: `figures/chapter-theory-of-sparse-regularization/certif-cs`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
