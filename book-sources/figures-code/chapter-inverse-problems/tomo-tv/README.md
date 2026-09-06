# Figure 9.10

Inverse Problems.

## Figure 9.10

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
TV-regularized reconstruction and the pseudoinverse for the same partial Fourier measurements.\relax
```

Minimizes .5||mask F x-y||²+.02 TV(x) for the same 13-line Cartesian Fourier acquisition as the pseudoinverse. The primal-dual method uses an exact Fourier-domain proximal step for fidelity and a pointwise Euclidean dual-ball projection for isotropic TV. The objective is checked against the initial pseudoinverse.

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--tomo-tv`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id inverse-problems--tomo-tv
```

Matching asset directory: `figures/chapter-inverse-problems/tomo-tv`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
