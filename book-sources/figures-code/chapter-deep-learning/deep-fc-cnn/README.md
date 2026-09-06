# Figure 17.1

Deep Learning.

## Figure 17.1

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Panel 1: fully connected network. Panel 2: convolutional neural network.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Retains the original three inputs, two hidden layers of four units, and three outputs. Every adjacent-layer pair is connected; line crossings are not nodes. The layer vectors are x0, x1, x2, x3, and component labels follow that notation. The output x3 denotes scores, avoiding confusion with the chapter target label y; a separate logistic map converts scores into probabilities. The dimension labels use the chapter notation: n_l=bar n_l*d_l and RGB input d_0=3. The five tensor blocks show an input, convolutional features, their downsampled version, a second feature stack, and a second downsampling. The activation arrows use the pointwise tilde-rho symbols because the chapter rho already includes downsampling. The output is the score vector x_L. Depth illustrates channel count and face size illustrates spatial resolution; exact numerical widths or filter sizes cannot be recovered and are not invented.

Omitted from the current comparison PDF. Stable identifier: `deep-learning--deep-fc-cnn`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id deep-learning--deep-fc-cnn
```

Matching asset directory: `figures/chapter-deep-learning/deep-fc-cnn`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
