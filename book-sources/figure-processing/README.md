# Figure comparisons

Open [figure-comparisons.pdf](figure-comparisons.pdf) to review all 105 original/TikZ
pairs. Each page contains:

- The exact figure number and complete caption from the main book.
- A link to the figure's page in `../FundationsDataScience.pdf`.
- A panel identifier when several reconstructions belong to one book figure.
- The original image, its reconstruction, and the mathematical interpretation notes.
- A stable drawing identifier in the footer.

Use the figure number and panel, or the stable identifier, when requesting edits.
The PDF bookmarks follow the book's figure order. Panel numbers refer to image
components in source order; the caption and reconstruction title identify their
meaning. These identifiers also appear above grouped panels in the main book.
Long groups may continue under the same figure number on another page; the
comparison link targets the page containing that particular panel.

All 105 reconstructions are now used in the main book and standalone chapters.
This separate volume preserves the original scans for comparison.

Run `make` from `book-sources/` to rebuild the book, all standalone chapters, and
this comparison PDF. `make book` rebuilds the book and comparison PDF together.
The book must be compiled first: the comparison captions are imported from its
LaTeX reference data, so they are not maintained as a separate copy.

`figure-index.json` records the current figure numbers, book pages, captions in
LaTeX, panel descriptions, stable IDs, asset paths, and publication status. The editable comparison
layout is in `comparison-layout.tex`; `figure-comparisons.tex` is its document
source. Drawing sources and interpretation notes remain in `../figures/tikz/`.
Build intermediates stay in `../build/figure-processing/`.
