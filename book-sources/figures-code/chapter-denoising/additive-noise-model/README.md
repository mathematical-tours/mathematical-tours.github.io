# Figure 7.2

Denoising.

## Figure 7.2

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Image acquisition, noise, and reconstruction by denoising.\relax
```

The clean signal and additive noise meet at a sum node, producing the observation Y. Only Y enters the denoiser D; the clean target remains unavailable to the estimator. The reconstruction is labeled as an estimate, not an exact inverse of the noisy acquisition.

Omitted from the current comparison PDF. Stable identifier: `denoising--additive-noise-model`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--additive-noise-model
```

Matching asset directory: `figures/chapter-denoising/additive-noise-model`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
