# Unnumbered illustration: nested spaces 1d

Wavelets.

## Unnumbered illustration: nested spaces 1d

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Unnumbered illustration: embeded spaces 1d
```

Restores the original two-row cascade: approximation spaces V_j on top and orthogonal detail spaces W_j below. Horizontal arrows are coarse projections; each diagonal extracts the detail complement from the preceding finer space. The identity underneath fixes the chapter convention that larger j is coarser. In two dimensions W_j^(2) contains all three tensor detail orientations.

Omitted from the current comparison PDF. Stable identifier: `wavelets--embeded-spaces-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--embeded-spaces-1d
```

Matching asset directory: `figures/chapter-wavelets/nested-spaces-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
