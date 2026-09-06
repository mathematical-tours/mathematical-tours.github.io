# Figure 5.21

Linear and Nonlinear Approximation.

## Figure 5.21

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
The same mesh over a lightened $P_1$ interpolant (left), with the boxed boundary region enlarged (right). Thin triangles follow the red singularity curve; the arrow indicates its tangent. Smooth regions use nearly equilateral elements.\relax
```

The same mesh as Figure 5.20 is shown over its smoothly varying nodal interpolant. The orange box is magnified with equal horizontal and vertical units. The long triangle edges track the red singularity curve; their short altitudes cross it. The orange arrow indicates the tangent direction. Away from the contour, the mesh is nearly equilateral, without radial spokes.

Omitted from the current comparison PDF. Stable identifier: `approximation--triangulation-anisoiso`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id approximation--triangulation-anisoiso
```

Matching asset directory: `figures/chapter-linear-and-nonlinear-approximation/triangulation-anisoiso`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
