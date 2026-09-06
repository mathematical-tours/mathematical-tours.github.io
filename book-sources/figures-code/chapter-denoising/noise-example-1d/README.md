# Figure 7.3

Denoising.

## Figure 7.3

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
One-dimensional additive Gaussian noise at low, medium, and high levels, with a common clean signal and normalized noise realization.\relax
```

Three noise levels use the same clean signal and the same standard Gaussian realization, scaled by sigma. This isolates noise amplitude from changes in the underlying sample. All panels share their vertical range.

Omitted from the current comparison PDF. Stable identifier: `denoising--noise-example-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--noise-example-1d
```

Matching asset directory: `figures/chapter-denoising/noise-example-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
