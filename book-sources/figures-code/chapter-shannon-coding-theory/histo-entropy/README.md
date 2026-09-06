# Figure 3.3

Shannon Coding Theory.

## Figure 3.3

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Three examples of probability distributions with corresponding entropies in bits.\relax
```

Recomputed the three normalized distributions and their base-two entropies. The uniform four-symbol distribution has exactly two bits of entropy. All bars share a vertical probability scale; the third entropy is rounded from its computed value.

Omitted from the current comparison PDF. Stable identifier: `shannon-coding--histo-entropy`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id shannon-coding--histo-entropy
```

Matching asset directory: `figures/chapter-shannon-coding-theory/histo-entropy`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
