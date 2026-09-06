# Figure 9.9

Inverse Problems.

## Figure 9.9

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Pseudoinverse reconstruction from partial Cartesian Fourier measurements on 13 and 32 radial lines.\relax
```

Implements the chapter’s partial-Fourier pseudoinverse model by setting unobserved Fourier entries to zero. The mask includes conjugate frequencies, so reconstruction is real. This is the exact pseudoinverse of that Cartesian restriction operator; it is not mislabeled filtered backprojection or the exact pseudoinverse of an interpolation-based discrete Radon matrix.

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--tomo-pinv`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id inverse-problems--tomo-pinv
```

Matching asset directory: `figures/chapter-inverse-problems/tomo-pinv`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
