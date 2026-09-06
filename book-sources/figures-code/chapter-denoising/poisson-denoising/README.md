# Figure 7.27

Denoising.

## Figure 7.27

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.28**.

Exact current book caption (LaTeX):

```tex
Poisson-corrupted flower images at increasing photon-count scales $\lambda _{\max }$.\relax
```

Replaces unspecified legacy count levels with a fully specified experiment. Photon counts are independent Poisson variables with mean lambda_max times the normalized Flower image; displayed counts are divided by lambda_max. Values above one are clipped only for display. The four acquisition parameters are now known and reproducible.

Omitted from the current comparison PDF. Stable identifier: `denoising--poisson-denoising`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--poisson-denoising
```

Matching asset directory: `figures/chapter-denoising/poisson-denoising`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
