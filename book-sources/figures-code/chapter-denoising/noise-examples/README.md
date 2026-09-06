# Figure 7.1

Denoising.

## Figure 7.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Gaussian, Poisson, and multiplicative Gamma noise simulated on the same flower image.\relax
```

These are controlled simulations of three acquisition-noise models on the same Flower photograph, not newly measured camera, microscopy, or radar data. The Gaussian standard deviation is 0.08, Poisson means are 30 times the normalized intensity, and the mean-one Gamma multiplier has shape 2. Display clipping does not modify the simulated arrays.

Omitted from the current comparison PDF. Stable identifier: `denoising--noise-examples`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--noise-examples
```

Matching asset directory: `figures/chapter-denoising/noise-examples`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
