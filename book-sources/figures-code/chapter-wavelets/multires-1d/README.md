# Figure 4.1

Wavelets.

## Figure 4.1

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Translations and dilations of a scaling function generate the approximation spaces. The Shannon scaling function $\phi (t)=\sin (\pi t)/(\pi t)$ is illustrated, with horizontal coordinate $x/2^j$ and heights relative to $2^{-j/2}$.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The chapter scaling convention is phi_(j,n)(x)=2^(-j/2) phi(2^(-j)x-n). The actual Shannon scaling function sinc is used, with common horizontal coordinate x/2^j and amplitude measured relative to 2^(-j/2). The three centers are 0, 4 times 2^j, and 0; their relative scales are 1, 2, and 1/2 and peaks 1, 1/sqrt(2), and sqrt(2). The third atom has translation index zero, as its caption requires.

Omitted from the current comparison PDF. Stable identifier: `wavelets--multires-1d`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--multires-1d
```

Matching asset directory: `figures/chapter-wavelets/multires-1d`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
