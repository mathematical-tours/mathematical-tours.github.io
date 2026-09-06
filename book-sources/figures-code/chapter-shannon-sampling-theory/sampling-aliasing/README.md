# Figure 1.4

Shannon Sampling Theory.

## Figure 1.4

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Sampling and Shannon interpolation in the frequency and time domains, with unit sample spacing. Left: the bandlimiting condition of Theorem~\ref {thm-shannon-sampling} holds. Right: overlapping spectral replicas cause aliasing; the reconstructed signal still passes through every sample. The rows show the original signal, sampling, the interpolation kernel, and reconstruction.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The paired layouts follow the supplied shannon-interp reference images. Unit sample spacing is used. Exact signal/spectrum pairs are sinc(t/a)^4 and a B3(a omega/(2 pi)), for a=4 and a=2. The rows show the original signal, its weighted Dirac comb and replicated spectrum, the sinc kernel and ideal low-pass filter, and the reconstruction. In the aliased case the bottom time curve is the exact inverse transform of the clipped replicated spectrum; it differs between samples while retaining every sample value. Red shading tracks folded spectral mass; the dashed red curves show the originals in the reconstruction row. The same vertical scales are used within each domain.

Omitted from the current comparison PDF. Stable identifier: `shannon--sampling-aliasing`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id shannon--sampling-aliasing
```

Matching asset directory: `figures/chapter-shannon-sampling-theory/sampling-aliasing`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
