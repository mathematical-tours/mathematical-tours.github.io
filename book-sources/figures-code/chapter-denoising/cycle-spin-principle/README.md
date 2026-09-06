# Unnumbered illustration: cycle spin principle

Denoising.

## Unnumbered illustration: cycle spin principle

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Unnumbered illustration: cycle spin principle
```

A LaTeX commutative square expresses translation equivariance: D T_tau = T_tau D, with T_tau f(x)=f(x-tau). Both paths act on the same f and use the same translation. This is the chapter property restored by cycle spinning.

Omitted from the current comparison PDF. Stable identifier: `denoising--cycle-spin-principle`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id denoising--cycle-spin-principle
```

Matching asset directory: `figures/chapter-denoising/cycle-spin-principle`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
