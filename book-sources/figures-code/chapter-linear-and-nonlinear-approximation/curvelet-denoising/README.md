# Figure 5.26

Linear and Nonlinear Approximation.

## Figure 5.26

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Comparison of translation-invariant wavelet denoising and curvelet denoising.\relax
```

Retains the original close-up/noisy/wavelet/curvelet comparison using a detail of the supplied flower. Both denoisers are recomputed from one seeded noisy observation, with oracle thresholds selected on the same grid and SNR measured against the known clean flower detail. Curvelets use the documented periodic oversampled Parseval frame, with smooth parabolic frequency windows and all pixel translations; this is not the fast wrapping implementation.

Omitted from the current comparison PDF. Stable identifier: `approximation--curvelet-denoising`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--curvelet-denoising
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/curvelet-denoising`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
