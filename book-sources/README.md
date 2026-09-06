# Mathematical Foundations of Data Sciences

The active `\include{chapters/...}` lines in `FundationsDataScience.tex` define the book and its 19 standalone chapters. Commented-out chapters and archived material are not included.

## Build

Install a full TeX Live distribution with `pdflatex` and `bibtex`, Python 3.12 or later, and the figure-generation dependencies:

```sh
python3 -m venv build/figure-runtime
build/figure-runtime/bin/python -m pip install -r figures-code/requirements.txt
make                    # book, all chapters, and complete figure audit
make book               # book and figure audit
python3 scripts/build_book.py --chapter fourier
```

The builder uses the isolated figure runtime automatically, or the interpreter specified by `FIGURE_PYTHON`. Python figure dependencies are pinned in `figures-code/requirements.txt`.

Outputs are `FundationsDataScience.pdf`, `chapters-pdf/<chapter>.pdf`, and `figure-processing/figure-comparisons.pdf`. Standalone chapters retain the book's chapter/equation/theorem numbering, contain their own references, and link to other chapters through `../FundationsDataScience.pdf`. Keep this relative layout when distributing the files.

The build works in `build/`, checks all final LaTeX logs, and stabilizes references before publishing. Warnings, missing characters and overfull/underfull boxes block publication. `build/build-report.json` records the result. Development previews may use `--allow-warnings`.

After an editorial change, run the full `make` before distributing the PDFs together. A selected-chapter build leaves other independent PDFs unchanged.

## Figures and author audit

All active figures have reproducible Python and/or TikZ sources in:

```text
figures/chapter-<chapter-title>/<figure-name>/
figures-code/chapter-<chapter-title>/<figure-name>/
```

**All 145 reviewed reconstructions are integrated into the book and all 19 independent chapters.** The 81 previously accepted TikZ displays are retained, along with the original cover, Figures 11.3 and 11.7, and Figures 15.8–15.15. Old Figure 6.6 is removed; old Figures 7.9 and 7.10 share one figure.

No figure comparisons remain pending. The [review status PDF](figure-processing/figure-comparisons.pdf) records completion; the [previous comparison volume](figure-processing/archive/2026-09-06-reviewed-figures/figure-comparisons.pdf) is archived with its matching book snapshot and numbering index. The complete inventory retains 239 entries in 236 figure folders, including the removed and merged entries.

Figures 5.20–5.22 use smoothly varying image regions and triangles elongated along their boundary. The mesh is nearly equilateral away from the edge, and its tangential alignment is checked numerically. Figure 8.2 uses stronger red/blue contrast while preserving the numerical gradient. Photographic inputs remain the supplied flower, the explicitly requested mandrill, and the documented 1919 Felix the Cat frame. New reconstructions do not use Lena.

See [generation instructions](figures-code/README.md), [audit instructions](figure-processing/README.md), and the [figure index](figure-processing/figure-index.md). `figure-processing/figure-project.json` is the detailed inventory; `figure-index.json` maps figure numbers to comparison pages and source directories. Source data are in `data/`; numerical provenance, seeds and checks accompany generated assets. Run `make figures` to regenerate figure assets and `make figure-audit` to verify the published project after a full build.

## Editorial conventions

All active chapter sources use UTF-8. The shared design is in `book-preamble.tex`, mathematical commands in `mystyle.sty`, and publication details in `book-metadata.tex`. Dates appear on the book cover and only on standalone chapter openings. The book-title/author/affiliation mastheads remain on independent chapter PDFs.

See the [final integration record](figure-processing/final-integration.md) and per-figure source guides for the figure corrections, acceptance decisions and validation results.
