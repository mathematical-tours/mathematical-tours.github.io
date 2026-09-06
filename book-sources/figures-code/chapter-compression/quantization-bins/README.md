# Unnumbered illustration: quantization bins

Compression.

## Unnumbered illustration: quantization bins

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Unnumbered illustration: quantization bins
```

Unnumbered dead-zone quantization diagram, in units of T. The zero cell is (-T,T). At |u|=T the symbol is nonzero; nonzero cells are reconstructed at their midpoints sign(q)(|q|+1/2)T, while zero is reconstructed exactly as zero. Endpoint dots distinguish the equality convention from strict hard thresholding.

Omitted from the current comparison PDF. Stable identifier: `compression--quantization-bins`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id compression--quantization-bins
```

Matching asset directory: `figures/chapter-compression/quantization-bins`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
