# Figure 7.13

Denoising.

## Figure 7.13

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.14**.

Exact current book caption (LaTeX):

```tex
Comparison of hard (left) and soft (right) thresholding.\relax
```

Shows the actual reconstructions at the hard/soft SNR maxima from the preceding plot, with each method tuned on the same grid and clean reference. Both displays use identical intensity bounds.

Omitted from the current comparison PDF. Stable identifier: `denoising--wavthresh-2d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--wavthresh-2d
```

Matching asset directory: `figures/chapter-denoising/wavthresh-2d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
