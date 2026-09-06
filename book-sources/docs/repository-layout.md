# Repository organization

The book is maintained in `book-sources/`, a subdirectory of the Mathematical Tours website repository. Run the commands below from that directory.

## Directory responsibilities

| Path | Responsibility |
| --- | --- |
| `FundationsDataScience.tex` | Main publication driver and authoritative chapter order. |
| `FundationsDataScience.pdf` | Published complete book; its historical filename is retained for existing links. |
| `chapters-tex/` | Active chapter text, included sections, shared configuration and bibliography. |
| `chapters-pdf/` | Published independent chapters only. |
| `figures/` | Chapter-organized figure assets, original panels and numerical provenance. |
| `figures-code/` | Matching Python/TikZ generators and pinned dependencies. |
| `data/` | Input images, audio and datasets with source records. |
| `tools/` | Maintained build, regeneration and audit utilities. |
| `docs/figures/` | Review PDFs, manifests, decisions, correction records and historical archives. |
| `aux/` | Auxiliary course notes, experiments and inactive drafts. |
| `legacy/` | Original MATLAB code, superseded text and historical migration utilities. |
| `build/` | Ignored local build products, caches, runtime and review evidence. |

## Cleanup on 6 September 2026

| Previous path | Current path |
| --- | --- |
| `chapters/` | `chapters-tex/` for active text; inactive drafts in `aux/chapters/`. |
| `book-metadata.tex`, `book-preamble.tex`, root `.sty` files | `chapters-tex/shared/`. |
| `biblio/` | `chapters-tex/bibliography/`. |
| `matlab/` | `legacy/matlab/`. |
| `auto-diff/`, `mle/`, `diffusion-models/`, `optim-ml/`, `writing-notes/` | Matching directories under `aux/`. |
| `scripts/` | `tools/`; superseded migration/correction utilities in `legacy/tools/`. |
| `figure-processing/` | `docs/figures/`. |
| `chapters/old/`, `chapters/perceptrons-old.tex` | `legacy/chapters/`. |
| `chapters/to-integrate/` | `aux/source-material/coding-blocks/`; content already incorporated in the book. |
| `chapters-pdf/compile_pdf.m` | `legacy/matlab/compile_pdf.m`. |
| `optimal-transport-OLD/` | Removed as requested, including deletion of its former tracked `optimal-transport/` paths. |

Tracked LaTeX intermediates and `.DS_Store` files were removed within `book-sources/`. Ignore rules now keep these files, Python caches and notebook checkpoints out of commits. Published PDFs, bibliography databases, input data and historical figure assets remain versioned. Already-deleted obsolete root `todo.tex` and `corrections.md` are not recreated; the figure correction record is [final-integration.md](figures/final-integration.md).

All active TeX inputs, figure assembly sources, manifests, build utilities and documentation use the new layout. Duplicate auxiliary copies of `mcode.sty` now use the shared package, whose declaration supports loading from a subdirectory. Historical archives retain the paths recorded when they were created. Archived review PDFs remain beside their matching historical book, while new reviews link to `../../FundationsDataScience.pdf` and standalone chapters continue to link to `../FundationsDataScience.pdf`.

## Maintenance

Use `make` to rebuild and audit the complete publication, `make check` to repeat its audit, and `make clean` to remove reproducible caches. Cleaning preserves the published PDFs, the local Python runtime and review evidence. Fresh builds seed figure references from checked-in assets, so no preexisting LaTeX auxiliary files are required.

Keep new active text in `chapters-tex/`, inputs in `data/`, matching figure assets/code in their chapter folders, and editorial records in `docs/`. Add draft material under `aux/`. Avoid putting scripts or temporary build outputs in `chapters-pdf/`.

## Validation

The build was run after clearing all active TeX and figure caches. The 277-page book, all 19 standalone chapters and the one-page review completion PDF compile with zero diagnostics. The publication audit passes all 767 shared-label checks, preserves the 515 historical asset hashes, and validates fonts, dates, figure acceptance and PDF destinations.

All 565 pages across the book and independent chapters were compared with their pre-cleanup PDFs: extracted text and page counts are unchanged. Raster comparison at 72 dpi found only five one-level antialiasing differences across three pages; visual inspection confirmed no content or layout changes. All 41 links in the latest archived figure review still resolve to its matching historical book.

The latest detailed record is [figures/verification.json](figures/verification.json).
