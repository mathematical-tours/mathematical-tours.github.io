# Figure 3.1

Shannon Coding Theory.

## Figure 3.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Left: the complete tree of binary words of length 3. Right: a prefix code.\relax
```

The full tree contains all eight length-three words. The pruned tree uses exactly the original code 0→001, 1→01, 2→1, 3→000. Edge labels append bits; the leaves show the source symbol and its codeword. The prefix property is checked.

Omitted from the current comparison PDF. Stable identifier: `shannon-coding--tree-prefix`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id shannon-coding--tree-prefix
```

Matching asset directory: `figures/chapter-shannon-coding-theory/tree-prefix`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
