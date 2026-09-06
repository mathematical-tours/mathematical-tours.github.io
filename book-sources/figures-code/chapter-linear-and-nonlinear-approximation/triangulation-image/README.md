# Figure 5.20

Linear and Nonlinear Approximation.

## Figure 5.20

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
A cartoon with smooth variation on both sides of an elliptical jump, a triangulation elongated along that boundary, and the continuous $P_1$ interpolant.\relax
```

The cartoon has different smooth nonaffine intensity functions inside and outside an elliptical jump. Nearly equilateral triangles fill the smooth regions; close contour-following layers create thin triangles whose long edges follow the tangent. Explicit strip connectivity prevents edge flips across the contour. The longest edges of all 128 jump-crossing triangles deviate by less than 5.4 degrees from the local tangent, with aspect ratios above 13.5. The third panel is the actual continuous P1 nodal interpolant. The contour is known analytically in this controlled example.

Omitted from the current comparison PDF. Stable identifier: `approximation--triangulation-image`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id approximation--triangulation-image
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/triangulation-image`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
