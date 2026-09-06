# Figure 12.6

Compressed Sensing.

## Figure 12.6

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Evolution of lower bounds $\hat \delta _k$ on the RIP constant.\relax
```

For each support size, the largest observed ||A_IᵀA_I-I||2 is computed over 160 sampled subsets. Cumulative maxima remain rigorous lower bounds by monotonicity of the true RIP constant. These sampled bounds cannot certify a small RIP constant and are not presented as exact values.

Omitted from the current comparison PDF. Stable identifier: `compressed-sensing--rip-const`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id compressed-sensing--rip-const
```

Matching asset directory: `figures/chapter-compressed-sensing/rip-const`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
