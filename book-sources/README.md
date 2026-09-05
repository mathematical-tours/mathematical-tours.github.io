# Mathematical Foundations of Data Sciences

The active `\include{chapters/...}` lines in `FundationsDataScience.tex` define
the book's contents and the list of standalone chapters. Commented-out chapters
and archived material are not part of the build.

## Build

Install Python 3 and a full TeX Live distribution (MacTeX on macOS), with
`pdflatex` and `bibtex` on your path. No Python packages are required.

```sh
make                 # complete book and every active standalone chapter
make book            # complete book only
python3 scripts/build_book.py --chapter fourier
```

The outputs are `FundationsDataScience.pdf` and `chapters-pdf/<chapter>.pdf`.
The selected-chapter command also rebuilds the book to obtain current reference
numbers. Each standalone chapter contains its own references, retains the book's
chapter/equation/theorem numbering, and links references to other chapters to
`../FundationsDataScience.pdf`. Keep that relative layout when distributing the
PDFs. Standalone page numbers start at 1.

After editing the book, run the full `make` before distributing the PDFs together.
The book-only and selected-chapter commands leave other chapter PDFs unchanged;
their reference numbers and link destinations may then be out of date.

The build runs in `build/`, without deleting or overwriting source-folder auxiliary
files. It reruns LaTeX until references stabilize and checks the final logs for
warnings, missing characters, and overfull or underfull boxes. A build with any
remaining diagnostics does not publish PDFs. The machine-readable report is
`build/build-report.json`; detailed logs are retained below `build/book/` and
`build/chapters/`. `--allow-warnings` is available for development previews.

All active chapter sources use UTF-8. The shared design is in `book-preamble.tex`;
mathematical commands are in `mystyle.sty`. See `corrections.md` for the editorial
and mathematical revision record.

The corrected classification-loss, Poisson variance-stabilization, and wavelet-support
figures have editable vector sources beside their PDFs. The build regenerates them
when their sources change; the Poisson data are computed by `scripts/variance_stabilization.py`.

The Gaussian-maximum standard-deviation plot uses a corrected PDF with its axis
label moved below the curve. The original empirical plot is preserved. Optional
regeneration uses `scripts/correct_max_gaussian_std.py` and requires `pypdf`;
the normal build uses the included PDF and needs no additional Python package.

## Reviewing the hand-drawn figures

The book and affected standalone chapters include **Figure reconstructions for
review** sections, accessible from the book's contents and PDF bookmarks. Each
landscape page places the original scan beside an editable TikZ reconstruction,
with a note explaining the mathematics, deliberate corrections, and any uncertain
interpretation. The original figures remain in the chapter text while the new
versions are reviewed.

There are 91 comparisons across 14 chapters. The sources, shared drawing style,
and inventory are in [`figures/tikz/`](figures/tikz/README.md). The normal build
compiles changed TikZ sources and regenerates the comparison sections before
building the PDFs. This uses TikZ, PGFPlots, and `standalone`, included in a full
TeX Live installation. The figure-by-figure correction record is in
[`corrections.md`](corrections.md#hand-drawn-figure-reconstruction-pass).
