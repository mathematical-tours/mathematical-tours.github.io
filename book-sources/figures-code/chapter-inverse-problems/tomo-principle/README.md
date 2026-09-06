# Figure 9.7

Inverse Problems.

## Figure 9.7

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
The modified Shepp--Logan phantom and one Radon projection at angle $\theta =\pi /6$. The marked line and point correspond to the same detector coordinate.\relax
```

The conventional modified Shepp-Logan test object is generated from ten ellipses. The red line has normal n_theta and signed offset t; the red dot is its actual line integral on the exact analytic one-dimensional Radon projection. Integrating the projection recovers the phantom mass. Ellipse conventions follow the standard modified Shepp-Logan model described in https://www.mathworks.com/help/images/ref/phantom.html .

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--tomo-principle`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id inverse-problems--tomo-principle
```

Matching asset directory: `figures/chapter-inverse-problems/tomo-principle`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
