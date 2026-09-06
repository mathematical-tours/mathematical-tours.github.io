# Figure 6.4

Compression.

## Figure 6.4

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
JPEG-2000 coding architecture.\relax
```

The architecture separates transform, quantization, embedded code-block coding, context-adaptive arithmetic coding, and rate/layer organization. It distinguishes irreversible CDF 9/7 coding from reversible integer 5/3 coding. The diagram describes the stages in the chapter without implying that arbitrary byte truncation remains decodable.

Omitted from the current comparison PDF. Stable identifier: `compression--jpeg-2k-architecture`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id compression--jpeg-2k-architecture
```

Matching asset directory: `figures/chapter-compression/jpeg-2k-architecture`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
