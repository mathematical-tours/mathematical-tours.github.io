# Figure 7.8

Denoising.

## Figure 7.8

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
SNR as a function of filter width for a signal (left) and an image (right).\relax
```

These curves are measured against the known clean signal and flower image for the same observations used in the neighboring examples. The marked grid maximizers are oracle choices for illustration; selecting them requires the clean reference and is not presented as a deployable estimator.

Omitted from the current comparison PDF. Stable identifier: `denoising--iltering-optimal-curve`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--iltering-optimal-curve
```

Matching asset directory: `figures/chapter-denoising/filtering-optimal-curves`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
