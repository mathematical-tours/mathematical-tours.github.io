# Editable reconstructions of hand-drawn figures

These 91 standalone TikZ drawings reconstruct the active book's hand-drawn
illustrations. The four JSON manifests map every source to its original asset,
owning chapter, mathematical context, and interpretation notes. The drawings use
the definitions and conventions of the current chapter text. Illustrative points
and shapes are identified as examples rather than measurements recovered from a
low-resolution scan.

Run `make` from the repository root to compile all drawings, the complete book,
and every standalone chapter. To compile one drawing directly, also from the
repository root:

```sh
pdflatex -output-directory=build figures/tikz/shannon-sinc.tex
```

Every `.tex` file is editable and includes `tikz-preamble.tex` for the shared
fonts, colors, and drawing styles. The adjacent PDFs are included in the book;
the build updates them when the drawing or its shared preamble changes.
Temporary figure logs are in `build/tikz-figures/`.

The generated `reviews/<chapter>.tex` files contain the side-by-side panels and
their notes. Edit the corresponding JSON manifest to change review prose;
`scripts/tikz_reviews.py` regenerates these files. Manifest prose is plain text,
not LaTeX. A `qualified` confidence value adds “Interpretation to check” to the
printed note. It identifies an unclear or schematic original, not a change to
the surrounding theorem's assumptions.
Book labels such as `eq-convol-1d` or `fig-fft` in the context are printed as
clickable equation or figure numbers, so readers can reach the supporting text.

The book includes each review section immediately after its chapter. The same
section is included when the chapter is compiled independently. Existing figure
numbers and original assets are preserved. The page footer carries the stable
manifest ID, so feedback can identify a reconstruction unambiguously.

The detailed inventory and mathematical decisions are recorded in
[`corrections.md`](../../corrections.md#hand-drawn-figure-reconstruction-pass).
