# Figure 1.9

Shannon Sampling Theory.

## Figure 1.9

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Quantizing an image with $K\in \{2,3,4,16\}$ gray levels. The normalized intensities lie in $[0,1)$; finite-level quantizers also specify clipping or endpoint conventions.\relax
```

Nearest-level quantization uses K equally spaced reproduction values between zero and one: Q(u)=clip(floor((K-1)u+1/2),0,K-1)/(K-1). All panels use the same flower image and display range. The finite endpoints and tie convention are explicit.

Omitted from the current comparison PDF. Stable identifier: `shannon--section3-quantize`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id shannon--section3-quantize
```

Matching asset directory: `figures/chapter-shannon-sampling-theory/section3-quantize`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
