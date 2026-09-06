# Figure 4.23

Wavelets.

## Figure 4.23

Accepted reconstruction, integrated into the book and independent chapters.

Previous audit identifier: **4.22**.

Exact current book caption (LaTeX):

```tex
One step of the 2-D wavelet transform algorithm.\relax
```

Restores the input photograph, one vertical split after row analysis, and the four quadrants after column analysis. LL/LH/HL/HH describe low/high filtering of the first and second array coordinates. One orthonormal Haar step preserves shape and energy. Coefficient displays use a shared signed logarithmic grayscale with zero at mid-gray; the coarse approximation is displayed at its image intensity scale. Colored boundaries follow the actual dyadic subband slices, from red at the finest level to blue, green and purple at coarser levels. Display compression does not alter the numerical coefficients.

Omitted from the current comparison PDF. Stable identifier: `wavelets--wavalgo-2d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--wavalgo-2d
```

Matching asset directory: `figures/chapter-wavelets/wavalgo-2d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
