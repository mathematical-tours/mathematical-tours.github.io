# Figure 13.15

Basics of Machine Learning.

## Figure 13.15

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Logistic classification in one and two dimensions, showing how $\norm {\be }$ controls the width of the probability transition.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The probability is theta(beta x)=1/(1+exp(-beta x)), with positive slopes 1 and 3. Both curves pass through probability 1/2 at x=0. Colored dots and guides locate the actual 0.1 and 0.9 crossings, and each dimension arrow measures the corresponding width 2 log(9)/beta. The larger slope has one third of the transition width. Both panels use beta=t v with the same unit normal v. The red dashed contours are labeled by their probabilities 0.1 and 0.9; the dimension arrows measure their normal separation, exactly 2 log(9)/t. The boundary at probability 1/2 stays fixed. Parameter scaling controls probability sharpness, not whether a data set is linearly separable.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--logistic-sharpness`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--logistic-sharpness
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/logistic-sharpness`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
