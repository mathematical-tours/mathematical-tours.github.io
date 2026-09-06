# Book text

The active `\include` list in `../FundationsDataScience.tex` defines chapter order. Files named `sec-*.tex` and `machine-learning-sec-pac.tex` are included sections, not additional standalone chapters.

- `shared/book-metadata.tex`: title, author, affiliation and edition date.
- `shared/book-preamble.tex`: fonts, layout, headers and shared LaTeX configuration.
- `shared/mystyle.sty`: mathematical commands; `mcode.sty` supplies code listings.
- `shared/notations_ot.sty`: additional notation retained for auxiliary material.
- `bibliography/`: BibTeX source databases.

Compile from `book-sources/` with `make`, or build one chapter with `python3 tools/build_book.py --chapter fourier`. TeX editor root directives point to the main driver. Build products belong in `build/`, and the published standalone PDFs belong in `chapters-pdf/`.

Inactive draft chapters are in `../aux/chapters/`; superseded text is in `../legacy/chapters/`.
