# Figure 1.8

Shannon Sampling Theory.

## Figure 1.8

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
A uniform scalar quantizer.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The axes are the normalized input u/T and integer index v. Filled left and open right endpoints implement the stated half-open cells exactly. The red vertical segment measures quantization error at a fixed input; rescaling gives the bound T/2 after dequantization. The diagonal represents the unquantized value.

Omitted from the current comparison PDF. Stable identifier: `shannon--shannon-quantizer`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id shannon--shannon-quantizer
```

Matching asset directory: `figures/chapter-shannon-sampling-theory/shannon-quantizer`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
