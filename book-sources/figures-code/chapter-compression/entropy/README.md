# Figure 6.3

Compression.

## Figure 6.3

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Three different probability distributions on the same alphabet of $Q=7$ symbols. Entropy is measured in bits.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. All three illustrations use the same seven-symbol alphabet (-3,...,3) and a common vertical scale. The probabilities are p_i=1/7 for every symbol. They are nonnegative and sum to one. Entropy is log_2(7), approximately 2.807 bits, using the convention 0 log(0)=0. These distributions are illustrative probability models, not empirical measurements. The probabilities are p_0=1 and p_i=0 otherwise. Entropy is zero bits, using the convention 0 log(0)=0. The probabilities are (0.04, 0.10, 0.20, 0.32, 0.20, 0.10, 0.04). Entropy is approximately 2.491 bits, using the convention 0 log(0)=0.

Omitted from the current comparison PDF. Stable identifier: `compression--entropy`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id compression--entropy
```

Matching asset directory: `figures/chapter-compression/entropy`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
