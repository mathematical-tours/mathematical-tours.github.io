# Figure 2.10

Fourier and Convolution.

## Figure 2.10

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Diagram of one radix-two FFT step. Both halves of the input enter the sum and difference; the outputs of the two smaller transforms are interleaved into even and odd frequencies.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Four distinct arrows connect both input halves to both the sum and difference. The difference is multiplied componentwise by the twiddle-factor vector before its half-size FFT. On the right, individual arrows map the two output streams into alternating even and odd positions; the first three entries are shown for N>=8. The zero-based indices match the DFT and the interleaving definition in the chapter. Crossings are not connections.

Omitted from the current comparison PDF. Stable identifier: `fourier--fft`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id fourier--fft
```

Matching asset directory: `figures/chapter-fourier-and-convolution/fft`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
