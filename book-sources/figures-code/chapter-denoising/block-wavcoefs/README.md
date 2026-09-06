# Figure 7.23

Denoising.

## Figure 7.23

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **7.24**.

Exact current book caption (LaTeX):

```tex
Left: wavelet coefficients, center: block thresholded coefficients, right: denoised image.\relax
```

Disjoint 4 by 4 blocks are formed separately within every orthogonal detail subband. Each block shares the Stein factor max(1-T^2/E_B,0), where E_B is its mean squared coefficient magnitude. The threshold is oracle-selected on the common grid, and the displayed coefficient arrays come from the same reconstruction. Coefficient displays use a shared signed logarithmic grayscale with zero at mid-gray; the coarse approximation is displayed at its image intensity scale. Colored boundaries follow the actual dyadic subband slices, from red at the finest level to blue, green and purple at coarser levels. Display compression does not alter the numerical coefficients.

Omitted from the current comparison PDF. Stable identifier: `denoising--block-wavcoefs`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--block-wavcoefs
```

Matching asset directory: `figures/chapter-denoising/block-wavcoefs`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
