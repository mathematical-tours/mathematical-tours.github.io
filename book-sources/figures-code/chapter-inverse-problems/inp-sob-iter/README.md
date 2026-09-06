# Figure 9.4

Inverse Problems.

## Figure 9.4

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Sobolev projected gradient descent algorithm.\relax
```

Every update is the actual projected Sobolev step x <- projection_C(x+.24 Delta x). Projection restores all known pixels exactly, and the step is below 1/4. Reported energies decrease and measured-data consistency is checked.

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--inp-sob-iter`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id inverse-problems--inp-sob-iter
```

Matching asset directory: `figures/chapter-inverse-problems/inp-sob-iter`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
