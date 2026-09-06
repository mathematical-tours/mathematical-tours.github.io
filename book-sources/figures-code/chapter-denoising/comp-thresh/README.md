# Figure 7.25

Denoising.

## Figure 7.25

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.26**.

Exact current book caption (LaTeX):

```tex
Left: translation invariant wavelet hard thresholding, center: block orthogonal Stein thresholding, right: block translation invariant Stein thresholding.\relax
```

All three methods are run and independently oracle-tuned on one observation. Stationary block Stein averages all 4 by 4 block alignments within each undecimated subband, so the combined estimator is equivariant to integer pixel translations. This property is checked numerically; averaging only a few image shifts is not labeled fully translation invariant.

Omitted from the current comparison PDF. Stable identifier: `denoising--comp-thresh`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--comp-thresh
```

Matching asset directory: `figures/chapter-denoising/comp-thresh`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
