# Figure 14.8

Optimization \& Machine Learning: Smooth Optimization.

## Figure 14.8

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Left: point clouds $(a_i)_i$ and their principal directions. Right: the quadratic part of $f(x)$.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Both panels share the eigenvectors u₁,u₂. The twelve reflection-symmetric points are centered and have empirical covariance eigenvalues λ₁=1 and λ₂=0.325 squared. Every covariance ellipse and quadratic level set is derived from these same eigenvalues, with reciprocal axis ratios. Width annotations sit above the panels so all data points and complete ellipses remain visible. The caption specifies the multiplier used for each family of semiaxes; C/n follows the covariance normalization in the text.

Omitted from the current comparison PDF. Stable identifier: `optim-ml-smooth--link-pca`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id optim-ml-smooth--link-pca
```

Matching asset directory: `figures/chapter-optimization-and-machine-learning-smooth-optimization/link-pca`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
