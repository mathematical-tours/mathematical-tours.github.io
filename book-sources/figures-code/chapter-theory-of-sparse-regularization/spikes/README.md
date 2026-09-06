# Figure 11.10

Theory of Sparse Regularization.

## Figure 11.10

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Convolution operator.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Two translated Gaussian kernels have unit peak; the spike heights and the output use the same vertical scale, so the displayed output is exactly their weighted sum. Kernel centers, spike locations, and output-location labels agree across the panels. The diagram distinguishes index-space weights from physical reconstruction-grid locations, and writes the output as both the continuous operator applied to the discrete measure and the dictionary product Ax.

Omitted from the current comparison PDF. Stable identifier: `sparse-theory--spikes`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id sparse-theory--spikes
```

Matching asset directory: `figures/chapter-theory-of-sparse-regularization/spikes`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
