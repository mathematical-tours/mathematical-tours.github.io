# Figure 1.3

Shannon Sampling Theory.

## Figure 1.3

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Spatial and temporal detail: a flower image with a pixelated crop, and a recorded birdsong waveform with a short-time enlargement. Matching red boxes identify the spatial crop.\relax
```

The red spatial box selects precisely the flower crop shown as 20 by 20 area-averaged pixels. The sound panels use the supplied bird recording and a 2 ms excerpt around a strong chirp; the red interval marks this exact excerpt. No temporal subsampling is applied.

Omitted from the current comparison PDF. Stable identifier: `shannon--discretization`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id shannon--discretization
```

Matching asset directory: `figures/chapter-shannon-sampling-theory/discretization`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
