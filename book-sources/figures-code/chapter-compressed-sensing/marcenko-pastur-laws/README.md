# Figure 12.3

Compressed Sensing.

## Figure 12.3

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Marchenko--Pastur densities $f_\be $ for several aspect ratios $\be $.\relax
```

Exact Marchenko–Pastur densities are evaluated only on their support [(1-sqrt(beta))²,(1+sqrt(beta))²]. All displayed beta values are below one, avoiding an unrepresented point mass at zero.

Omitted from the current comparison PDF. Stable identifier: `compressed-sensing--marcenko-pastur-laws`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id compressed-sensing--marcenko-pastur-laws
```

Matching asset directory: `figures/chapter-compressed-sensing/marcenko-pastur-laws`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
