# Figure 12.1

Compressed Sensing.

## Figure 12.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Simulated compressed sensing of the flower image, with $Q=4096$ pixels. The three wavelet-$\ell ^1$ reconstructions use $P=512$, $1024$, and $2048$ orthonormal signed Walsh--Hadamard measurements, respectively.\relax
```

Simulates compressed sensing of the supplied flower with three nested measurement counts. The sensing rows are orthonormal Walsh-Hadamard patterns with a fixed random pixel permutation; signed measurements can be formed by complementary exposures. Every reconstruction solves the stated wavelet-synthesis l1 problem from its own measurements using the same lambda and iterations. The operator and adjoint are exact, with norm one. These are simulated acquisitions, not new camera measurements.

Omitted from the current comparison PDF. Stable identifier: `compressed-sensing--single-pixel`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id compressed-sensing--single-pixel
```

Matching asset directory: `figures/chapter-compressed-sensing/single-pixel`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
