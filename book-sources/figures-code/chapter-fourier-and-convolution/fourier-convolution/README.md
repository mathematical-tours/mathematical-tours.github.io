# Figure 2.3

Fourier and Convolution.

## Figure 2.3

Previously accepted TikZ drawing retained in the reading editions.

Exact current book caption (LaTeX):

```tex
Convolution on $\RR $. The dashed curve is $f\star g$; the marked value $(f\star g)(x)$ equals the shaded integral.\relax
```

Regenerated from the accepted editable TikZ figure; its geometry is retained. The shaded interval is [x-epsilon,x+epsilon], the support of g(x-t), with g the unnormalized box on [-epsilon,epsilon]. The dashed red curve is the exact sliding integral of the same displayed signal, computed by integrating its affine and sinusoidal terms; the red dot is at (x,(f*g)(x)). The illustrated width is 2 epsilon=1.5, so the convolution is an integral, not a unit-mass average. The curve and shaded area use the same normalization.

Omitted from the current comparison PDF. Stable identifier: `fourier--fourier-convolution`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python tools/regenerate_figures.py --id fourier--fourier-convolution
```

Matching asset directory: `figures/chapter-fourier-and-convolution/fourier-convolution`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

The TikZ sources and `proposed.tex` specify the editable drawing and its assembly.
