# Figure 1.2

Shannon Sampling Theory.

## Figure 1.2

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
A color flower image and a synthetic multispectral image with 32 channels, with an example spectrum.\relax
```

The RGB panel is the supplied flower photograph. The 32-band cube is a schematic synthetic spectral field with two wavelength-dependent components; it is not inferred from RGB channels and does not claim real multispectral measurements. This separates the spatial dimension d=2 from the channel count s.

Omitted from the current comparison PDF. Stable identifier: `shannon--examples-2`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id shannon--examples-2
```

Matching asset directory: `figures/chapter-shannon-sampling-theory/examples-2`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
