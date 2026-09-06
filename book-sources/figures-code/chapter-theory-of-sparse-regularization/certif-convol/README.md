# Figure 11.11

Theory of Sparse Regularization.

## Figure 11.11

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Continuous representations of the precertificates $\eta _F$ and $\eta _V$ for a convolution operator $A$.\relax
```

Gaussian convolution gives C(x,z)=exp(-(x-z)²/(4 sigma²)), up to an irrelevant common normalization. Both precertificates are obtained from the kernel systems in the text, including the mixed derivative block with its correct sign. Derivative constraints enforce stationarity at support points; feasibility is measured separately and is not assumed.

Omitted from the current comparison PDF. Stable identifier: `sparse-theory--certif-convol`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id sparse-theory--certif-convol
```

Matching asset directory: `figures/chapter-theory-of-sparse-regularization/certif-convol`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
