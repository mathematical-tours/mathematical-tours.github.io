# Figure 7.20

Denoising.

## Figure 7.20

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.21**.

Exact current book caption (LaTeX):

```tex
Left: semi-soft thresholder, right: Stein thresholder.\relax
```

Plots the exact semi-soft rules with transition interval [T,mu T] and the exact Stein attenuation max(1-T^2/x^2,0)x, defining the origin to be zero. Soft thresholding is shown for comparison: its large-amplitude shrinkage tends to T, whereas Stein shrinkage tends to zero.

Omitted from the current comparison PDF. Stable identifier: `denoising--semisoft`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--semisoft
```

Matching asset directory: `figures/chapter-denoising/semisoft`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
