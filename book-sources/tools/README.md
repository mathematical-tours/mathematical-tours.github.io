# Publication tools

Run these commands from `book-sources/`. `make help` lists the common entry points.

| Tool | Purpose |
| --- | --- |
| `build_book.py` | Build the book, standalone chapters and current figure review PDF. |
| `regenerate_figures.py` | Rebuild figures from Python/TikZ sources and track dependencies. |
| `build_figure_comparator.py` | Extract exact book numbers/captions and build pending comparisons. |
| `document_figure_project.py` | Refresh per-figure guides, review index and author decisions. |
| `audit_figure_project.py` | Check accepted assets, data hashes, fonts, labels, dates and PDF links. |
| `figure_inventory.py` | Inventory figures in active chapters and their included sections. |
| `figure_tex.py` | Shared parsing of TeX groups and stabilized book references. |
| `clean_build.py` | Remove rebuildable caches while preserving PDFs and the Python runtime. |

The figure runtime is configured through `FIGURE_PYTHON`, defaulting to `build/figure-runtime/bin/python`. Dependencies are pinned in `../figures-code/requirements.txt`.

A fresh build first obtains reference labels from the checked-in reading assets, then regenerates figures and compiles the final editions. Figure caches include the shared metadata, preamble and style dependencies. `make` finishes with the publication audit and refreshes figure documentation; `make book` and selected-chapter builds produce partial publications and should be followed by a full `make` before distribution.

One-time migrations and superseded figure utilities are preserved in `../legacy/tools/` as historical reference. They are not part of the current pipeline.
