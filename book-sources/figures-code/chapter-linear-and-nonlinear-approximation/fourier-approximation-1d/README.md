# Figure 5.13

Linear and Nonlinear Approximation.

## Figure 5.13

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Fourier approximation of a signal.\relax
```

Low-frequency real Fourier projections of the same piecewise-smooth periodic signal show Gibbs oscillations near jumps. Counts include the constant and individual sine/cosine atoms; ties in frequency are resolved deterministically.

Omitted from the current comparison PDF. Stable identifier: `approximation--ourier-approx-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--ourier-approx-1d
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/fourier-approximation-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
