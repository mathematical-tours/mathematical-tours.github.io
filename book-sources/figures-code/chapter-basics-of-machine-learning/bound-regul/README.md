# Figure 13.11

Basics of Machine Learning.

## Figure 13.11

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Conditional expectation.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Here X is uniform on [-1.75,1.75] and Y=m(X)+epsilon, with independent centered asymmetric mixture noise: weights 3/4 and 1/4, Gaussian means -0.3 and 0.9, and standard deviation 0.28. Shading encodes joint density; dividing its vertical slice by p_X(x0) gives the right-hand conditional density. Its mean m(x0) differs from its mode. The theorem does not assume this illustrative noise model.

Omitted from the current comparison PDF. Stable identifier: `machine-learning--bound-regul`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id machine-learning--bound-regul
```

Matching asset directory: `figures/chapter-basics-of-machine-learning/bound-regul`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
