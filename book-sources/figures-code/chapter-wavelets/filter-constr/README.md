# Figure 4.25

Wavelets.

## Figure 4.25

Previously accepted TikZ drawing retained in the reading editions.

Previous audit identifier: **4.24**.

Exact current book caption (LaTeX):

```tex
Constraints on low-pass and high-pass filters.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. Squared magnitudes use the order-four Daubechies half-band polynomial: s=sin^2(omega/2), H2=2(1-s)^4(1+4s+10s^2+20s^3), G2=2-H2=H2(omega+pi). They are exactly complementary and have high-order flat extrema. The phases can be chosen by spectral factorization and the quadrature-mirror relation; magnitudes alone do not express the phase cancellation in C4. This merges the former duplicate complementarity figure into the filter-constraints figure.

Omitted from the current comparison PDF. Stable identifier: `wavelets--filter-constr`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id wavelets--filter-constr
```

Matching asset directory: `figures/chapter-wavelets/filter-constr`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
