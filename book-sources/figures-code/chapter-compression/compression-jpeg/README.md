# Figure 6.5

Compression.

## Figure 6.5

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
JPEG and JPEG-2000 reconstructions of the flower (first pair) and mandrill (second pair), at matched target bit rates.\relax
```

One row compares JPEG and JPEG-2000 on the flower and a different second image, the author-requested mandrill. Both encoders share the same uncompressed reference and a 0.5 bits/pixel target. Labels count actual stream bytes including headers; panels are decoded from the saved streams.

Omitted from the current comparison PDF. Stable identifier: `compression--compression-jpeg`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id compression--compression-jpeg
```

Matching asset directory: `figures/chapter-compression/compression-jpeg`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
