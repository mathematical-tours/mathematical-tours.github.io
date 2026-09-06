# Figure 4.18

Wavelets.

## Figure 4.18

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.17**.

Exact current book caption (LaTeX):

```tex
2-D wavelets and a schematic support centered at $(2^j n_1,2^j n_2)$, with width $K2^j$ (right).\relax
```

Restores the three surface plots and the separate support enclosure of the original. These are actual tensor products of the Daubechies-3 scaling and wavelet functions, with shared height and spatial scales. An integer shift puts the generator support in [-2,3], hence in the valid centered K=6 enclosure; the enclosure is not claimed to be minimal.

Omitted from the current comparison PDF. Stable identifier: `wavelets--wavelets-2d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id wavelets--wavelets-2d
```

Matching asset directory: `figures/chapter-wavelets/wavelets-2d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
