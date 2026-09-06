# Figure 9.5

Inverse Problems.

## Figure 9.5

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Sobolev inpainting of a grid-shaped occlusion on the flower image.\relax
```

Recreates the cage-removal setup using the supplied flower photograph and an explicit grid-shaped occlusion. Both the known-pixel mask and the clean reference are available for this controlled experiment. The Sobolev reconstruction enforces observed pixels exactly after every step; no legacy parrot image is used.

Omitted from the current comparison PDF. Stable identifier: `inverse-problems--inpainting-parrot`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id inverse-problems--inpainting-parrot
```

Matching asset directory: `figures/chapter-inverse-problems/inpainting-parrot`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
