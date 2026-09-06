# Figure 3.6

Shannon Coding Theory.

## Figure 3.6

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Successive merges in Huffman's coding algorithm.\relax
```

The merges are computed from probabilities (0.27,0.53,0.19,0.01), rather than inferred from the old drawing. Red bars identify the two least likely current nodes. The resulting tree gives code lengths (2,1,3,3), expected length 1.67 bits, and the correct entropy bound.

Omitted from the current comparison PDF. Stable identifier: `shannon-coding--huff-code`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id shannon-coding--huff-code
```

Matching asset directory: `figures/chapter-shannon-coding-theory/huff-code`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
