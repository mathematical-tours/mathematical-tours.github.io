# Figure 13.12

Basics of Machine Learning.

## Figure 13.12

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Linear regression.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The red line is the least-squares affine fit to the actual displayed illustrative points. The dashed curve is the model conditional mean used to construct those points. Residuals are vertical output errors, not orthogonal distances to the line.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--ml-linear-fit`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id machine-learning--ml-linear-fit
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/ml-linear-fit`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
