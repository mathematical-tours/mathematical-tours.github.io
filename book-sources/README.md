# Mathematical Foundations of Data Sciences

The active `\include{chapters-tex/...}` lines in `FundationsDataScience.tex` define the book and its 19 standalone chapters. Auxiliary notes and archived material are not included.

## Repository layout

| Directory | Contents |
| --- | --- |
| `chapters-tex/` | Active chapters, included sections, shared styles and bibliography. |
| `chapters-pdf/` | Published independent chapters. |
| `figures/` and `figures-code/` | Matching chapter-organized assets and Python/TikZ sources. |
| `data/` | Input images, audio and datasets. |
| `tools/` | Build, regeneration and audit utilities. |
| `docs/` | Editorial records and figure comparisons under `docs/figures/`. |
| `aux/` | Auxiliary course notes, experiments and inactive drafts. |
| `legacy/` | Original MATLAB code and superseded sources/tools. |
| `build/` | Ignored build products, caches and local runtime. |

The root keeps the main book driver and published PDF. See the [layout and migration guide](docs/repository-layout.md) for conventions and the cleanup record.

## Build

Run commands from this `book-sources/` directory. Install a full TeX Live distribution with `pdflatex` and `bibtex`, Python 3.12 or later, and the figure-generation dependencies:

```sh
python3 -m venv build/figure-runtime
build/figure-runtime/bin/python -m pip install -r figures-code/requirements.txt
make                    # book, all chapters, and complete figure audit
make book               # book and figure review PDF
make check              # audit a complete publication
make clean              # clear caches; keep published PDFs and Python runtime
python3 tools/build_book.py --chapter fourier
```

The builder uses the isolated figure runtime automatically, or the interpreter specified by `FIGURE_PYTHON`. Python figure dependencies are pinned in `figures-code/requirements.txt`.

Outputs are `FundationsDataScience.pdf`, `chapters-pdf/<chapter>.pdf`, and `docs/figures/figure-comparisons.pdf`. Standalone chapters retain the book's chapter/equation/theorem numbering, contain their own references, and link to other chapters through `../FundationsDataScience.pdf`. Keep this relative layout when distributing the files.

The build works in `build/`, checks all final LaTeX logs, and stabilizes references before publishing. Warnings, missing characters and overfull/underfull boxes block publication. `build/build-report.json` records the result. Development previews may use `--allow-warnings`.

After an editorial change, run the full `make` before distributing the PDFs together. A selected-chapter build leaves other independent PDFs unchanged.

## Figures and author audit

All active figures have reproducible Python and/or TikZ sources in:

```text
figures/chapter-<chapter-title>/<figure-name>/
figures-code/chapter-<chapter-title>/<figure-name>/
```

**All 145 reviewed reconstructions are integrated into the book and all 19 independent chapters.** The 81 previously accepted TikZ displays are retained, along with the original cover, Figures 11.3 and 11.7, and Figures 15.8–15.15. Old Figure 6.6 is removed; old Figures 7.9 and 7.10 share one figure.

No figure comparisons remain pending. The [review status PDF](docs/figures/figure-comparisons.pdf) records completion; the [previous comparison volume](docs/figures/archive/2026-09-06-reviewed-figures/figure-comparisons.pdf) is archived with its matching book snapshot and numbering index. The complete inventory retains 239 entries in 236 figure folders, including the removed and merged entries.

Figures 5.20–5.22 use smoothly varying image regions and triangles elongated along their boundary. The mesh is nearly equilateral away from the edge, and its tangential alignment is checked numerically. Figure 8.2 uses stronger red/blue contrast while preserving the numerical gradient. Photographic inputs remain the supplied flower, the explicitly requested mandrill, and the documented 1919 Felix the Cat frame. New reconstructions do not use Lena.

See [generation instructions](figures-code/README.md), [audit instructions](docs/figures/README.md), and the [figure index](docs/figures/figure-index.md). `docs/figures/figure-project.json` is the detailed inventory; `docs/figures/figure-index.json` maps figure numbers to comparison pages and source directories. Source data are in `data/`; numerical provenance, seeds and checks accompany generated assets. Run `make figures` to regenerate figure assets and `make figure-audit` to verify the published project after a full build.

## Editorial conventions

All active chapter sources use UTF-8. The shared design is in `chapters-tex/shared/book-preamble.tex`, mathematical commands in `chapters-tex/shared/mystyle.sty`, and publication details in `chapters-tex/shared/book-metadata.tex`. Dates appear on the book cover and only on standalone chapter openings. The book-title/author/affiliation mastheads remain on independent chapter PDFs.

See the [final integration record](docs/figures/final-integration.md) and per-figure source guides for the figure corrections, acceptance decisions and validation results.
