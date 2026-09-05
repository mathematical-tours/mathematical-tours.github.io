# Editable reconstructions of hand-drawn figures

These 91 standalone TikZ drawings reconstruct the active book's hand-drawn
illustrations. The four JSON manifests map every source to its original asset,
owning chapter, mathematical context, and interpretation notes. The drawings use
the definitions and conventions of the current chapter text. Illustrative points
and shapes are identified as examples rather than measurements recovered from a
low-resolution scan.

Run `make` from the repository root to compile all drawings, the complete book,
every standalone chapter, and the separate comparison PDF. To compile one drawing directly, also from the
repository root:

```sh
pdflatex -output-directory=build figures/tikz/shannon-sinc.tex
```

Every `.tex` file is editable and includes `tikz-preamble.tex` for the shared
fonts, colors, and drawing styles. The adjacent PDFs are included in the book,
standalone chapters, and comparison volume;
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

The generated pages are collected in `figure-processing/figure-comparisons.pdf`;
they are not included in the book or standalone chapters. Each manifest's
`book_label` points to the numbered figure (or its continued caption), and `book_panel` identifies
its component when appropriate. The comparison heading imports the exact number
and full caption from the main book's stabilized reference data. The four formerly
unnumbered sketches now have captions in the main text. The page footer retains
the stable manifest ID, so feedback can identify a reconstruction unambiguously
even if later book edits change its number. The machine-readable mapping is in
`figure-processing/figure-index.json`.

The detailed inventory and mathematical decisions are recorded in
[`corrections.md`](../../corrections.md#hand-drawn-figure-reconstruction-pass).

All 91 manifest entries have `in_book: true`. The reading editions use the TikZ
PDFs at their natural size, reduced only when needed to fit the text width.
Separate reconstructed components have explicit panel headings. Four long
sequences use continued figures to retain their number and legible lettering;
the comparison links target the page containing the corresponding component.
Original assets remain available for review and are not deleted.

The four scalar soft-thresholding drawings are also reused in the advanced
optimization chapter. Both occurrences use the same TikZ PDFs.
