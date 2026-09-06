# Figure 4.13

Wavelets.

## Figure 4.13

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Analysis and inverse wavelet transforms in block-matrix form. Synthesis uses the transpose of the orthogonal analysis operator.\relax
```

Recreates the original block-operator picture: analysis has two vertically stacked convolution/downsampling blocks; synthesis transposes this arrangement into two horizontally concatenated upsampling/convolution blocks. The stacked vector orders a_j above d_j, and both sides act on the correctly sized a_(j-1). The concrete eight-sample Haar matrix is retained as a separate accepted figure.

Omitted from the current comparison PDF. Stable identifier: `wavelets--inversion-transpose`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--inversion-transpose
```

Matching asset directory: `figures/chapter-wavelets/inversion-transpose`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
