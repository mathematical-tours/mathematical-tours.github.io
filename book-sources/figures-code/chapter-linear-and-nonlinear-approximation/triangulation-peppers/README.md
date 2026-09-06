# Figure 5.22

Linear and Nonlinear Approximation.

## Figure 5.22

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Compression of a piecewise-smooth cartoon using a contour-aligned mesh and JPEG-2000 at the same transmitted bit budget. From left to right: original, mesh, decoded $P_1$ image, and JPEG-2000 image. Mesh costs include coordinates, nodal values, connectivity, and headers.\relax
```

Uses smooth nonaffine intensities on both sides of an elliptical jump. A quasiuniform mesh in smooth regions transitions to thin, tangential triangles along the boundary; the same 1390-triangle mesh is used in Figures 5.20 and 5.21. The continuous P1 reconstruction is decoded from a real mesh stream containing a 10-byte header, uint16 coordinates, uint8 nodal values and uint16 triangle indices. The mesh and JPEG2000 streams have equal transmitted byte lengths, including all connectivity, headers and padding. No optimal entropy-coding claim is made.

Omitted from the current comparison PDF. Stable identifier: `approximation--triangulation-peppers`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--triangulation-peppers
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/triangulation-peppers`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
