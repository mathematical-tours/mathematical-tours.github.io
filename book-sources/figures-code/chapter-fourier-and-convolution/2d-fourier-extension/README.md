# Figure 2.14

Fourier and Convolution.

## Figure 2.14

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
2D Fourier orthogonal bases.\relax
```

The complex Fourier basis is indexed by integer lattice points. The image panels show exact real parts cos(2 pi (k1 x1+k2 x2)); their orientation and oscillation count follow the marked lattice vectors. All panels have the same spatial domain and amplitude range.

Omitted from the current comparison PDF. Stable identifier: `fourier--2d-fourier-extension`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id fourier--2d-fourier-extension
```

Matching asset directory: `figures/chapter-fourier-and-convolution/2d-fourier-extension`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
