# Figure 3.5

Shannon Coding Theory.

## Figure 3.5

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Panel 1: full binary tree obtained by completing the tree associated with the code $(c_1=0, c_2=10, c_3=110, c_4=111)$. Panel 2: subtrees with the prescribed codeword lengths packed from left to right.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The code is (0,10,110,111). Its codewords reserve 4, 2, 1, and 1 leaves in the full depth-three tree. Tinted subtrees connect the highlighted codeword roots to disjoint descendant blocks. Edge bits show how a path spells a codeword; the formula explains the count 2^(m-l). For prescribed lengths (1,2,3), the blocks contain 4, 2, and 1 leaves. The left-to-right packing arrow and tinted subtrees expose the construction: choose the largest block first, then continue at the next dyadic boundary. This gives codewords (0,10,110), leaves 111 unused, and demonstrates 7/8 <= 1 without assuming equality.

Omitted from the current comparison PDF. Stable identifier: `shannon-coding--kraft-trees`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id shannon-coding--kraft-trees
```

Matching asset directory: `figures/chapter-shannon-coding-theory/kraft-trees`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
