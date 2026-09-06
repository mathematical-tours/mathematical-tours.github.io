# Figure 3.7

Shannon Coding Theory.

## Figure 3.7

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Difference coding. From left to right: the signal, its first differences, the corresponding histograms, and a code tree for the differences.\relax
```

A fixed discrete signal is differenced, and both histograms and the Huffman tree are recomputed from its actual counts. Reconstruction retains the initial value and cumulatively sums the differences. Difference symbols are more concentrated at zero for this slowly varying example; this is an example, not a universal decrease of marginal entropy.

Omitted from the current comparison PDF. Stable identifier: `shannon-coding--codage-differences`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id shannon-coding--codage-differences
```

Matching asset directory: `figures/chapter-shannon-coding-theory/codage-differences`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
