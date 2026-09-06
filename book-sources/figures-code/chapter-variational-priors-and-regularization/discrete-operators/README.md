# Figure 8.2

Variational Priors and Regularization.

## Figure 8.2

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Discrete gradient of the flower image. Red and blue encode the signed horizontal and vertical components with enhanced display contrast; neutral gray represents zero. The last panel shows the Euclidean gradient magnitude.\relax
```

The middle RGB image encodes signed horizontal differences in red and signed vertical differences in blue. A shared contrast curve sign(z) min(|z|/s,1)^0.45, with s equal to the pooled 98.5th absolute-gradient percentile, strengthens weak components and clips only the display tails. Red and blue channels equal 0.5 plus half this signed curve; green stays at 0.5. The gradient arrays themselves are unchanged. Zero gradient is neutral gray; increasing/decreasing channels distinguish signs. The final panel is the actual Euclidean magnitude. Periodic div=-gradient-adjoint is verified.

Omitted from the current comparison PDF. Stable identifier: `variational-priors--discrete-operators`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id variational-priors--discrete-operators
```

Matching asset directory: `figures/chapter-variational-priors-and-regularization/discrete-operators`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
