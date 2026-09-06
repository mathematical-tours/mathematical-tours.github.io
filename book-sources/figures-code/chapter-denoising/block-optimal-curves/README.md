# Figure 7.24

Denoising.

## Figure 7.24

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.25**.

Exact current book caption (LaTeX):

```tex
SNR as a function of $T/\si $ (left) and comparison of different block sizes (right).\relax
```

The left panel applies hard, soft, and Stein factors to the block RMS, with the block-size normalization included. The right panel independently optimizes T/sigma for each block side length. Blocks remain inside one scale and orientation, and every comparison uses the same noisy image.

Omitted from the current comparison PDF. Stable identifier: `denoising--block-optimal-curves`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id denoising--block-optimal-curves
```

Matching asset directory: `figures/chapter-denoising/block-optimal-curves`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
