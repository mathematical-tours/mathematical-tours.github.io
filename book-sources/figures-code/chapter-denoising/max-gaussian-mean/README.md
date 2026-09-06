# Figure 7.15

Denoising.

## Figure 7.15

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.16**.

Exact current book caption (LaTeX):

```tex
Empirical mean (top) and standard deviation (bottom) of the largest absolute Gaussian noise coefficient.\relax
```

For each N, computes max_i |Z_i| over 384 independent standard-Gaussian vectors and reports its empirical mean and standard deviation. Prefixes share draws across N only to smooth the comparison; trials remain independent. The dashed leading-order scale is sqrt(2 log N). This is the absolute maximum, unlike the signed maximum in the compressed-sensing chapter.

Omitted from the current comparison PDF. Stable identifier: `denoising--max-gaussian-mean`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--max-gaussian-mean
```

Matching asset directory: `figures/chapter-denoising/max-gaussian-mean`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
