# Figure 1.1

Shannon Sampling Theory.

## Figure 1.1

Accepted reconstruction, integrated into the book and independent chapters.

Exact current book caption (LaTeX):

```tex
Examples of sound ($d=1$), an image ($d=2$), and a video ($d=3$).\relax
```

The waveform uses the supplied bird recording with its actual sample rate and duration. The grayscale image is the supplied flower image. The video is an explicitly simulated sequence of translated flower frames; it illustrates d=3, not an additional measured recording.

Omitted from the current comparison PDF. Stable identifier: `shannon--examples-1`.

Rebuild from the repository root:

```sh
build/figure-runtime/bin/python scripts/regenerate_figures.py --id shannon--examples-1
```

Matching asset directory: `figures/chapter-shannon-sampling-theory/examples-1`. `context.tex` records the mathematical context at reconstruction time. `original/` preserves historical assets and `original.pdf` assembles the previous figure. `proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.

`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.
