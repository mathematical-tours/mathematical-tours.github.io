# Historical figure material

`legacy-assets/` is the preserved pre-migration `figures/` tree. It contains all historical assets and source code, including material not used by the active book. The active chapter sources now point to byte-identical files in chapter/figure `original/` directories.

`previous-figure-comparisons.pdf`, its index, source and layout record the earlier 105-panel scan-to-TikZ project. Those TikZ drawings were accepted before the present regeneration project. The new audit compares the current reading edition with the next proposals.

The archived comparison PDF retains its historical relative book links. To follow those links, temporarily place a copy one directory above, in `figure-processing/`, beside the current audit. The published current audit uses the correct layout directly.

Historical files are preserved for provenance; they are not the active generation entry points. Use `scripts/build_book.py` and `scripts/regenerate_figures.py` for the current build.
