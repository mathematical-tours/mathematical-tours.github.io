# Figure 5.7

Linear and Nonlinear Approximation.

## Figure 5.7

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
A grayscale frame from \emph {Feline Follies} (1919), a mathematical jump along a smooth graph, and a Gaussian-smoothed cartoon. Film frame: Pat Sullivan, via Wikimedia Commons, public domain.\relax
```

The first panel is a real grayscale animation still from Feline Follies (1919), attributed by Wikimedia Commons to Pat Sullivan and listed there as public domain. Source: https://commons.wikimedia.org/wiki/File:Felix_1919.jpg . The other panels illustrate a mathematically defined jump along a smooth graph and an actual Gaussian-smoothed cartoon. An animation drawing and the mathematical cartoon class are distinguished.

Omitted from the current comparison PDF. Stable identifier: `approximation--sample-cartoon`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--sample-cartoon
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/sample-cartoon`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
