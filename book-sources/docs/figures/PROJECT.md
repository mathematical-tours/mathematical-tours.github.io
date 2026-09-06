# Figure regeneration and publication

All reviewed figures are now published. The 239 inventory entries comprise 145 accepted reconstructions, 81 retained TikZ displays, 11 explicitly retained originals, one removed display and one merged display. There are 237 active displays in 236 historical/current figure folders. Chapter 16 has no figures.

The main book has 277 pages; all 19 independent chapters and the one-page review-completion PDF compile without diagnostics. No comparisons remain pending. The previous 43-entry audit is archived with its matching book and valid historical links under `archive/2026-09-06-reviewed-figures/`.

## Reproducibility

Asset and code directories use matching `chapter-<title>/<figure-name>/` paths under `figures/` and `figures-code/`. Historical assets remain in `original/` and their composite `original.pdf`; all 515 historical reading assets are byte-preserved. Accepted reconstructions keep the stable name `proposed.pdf`. Numerical generators record input hashes, random seeds, parameters and mathematical checks.

`tools/regenerate_figures.py` rebuilds accepted and pending reconstructions, skipping explicitly retained originals and removed/merged entries. `tools/build_book.py` generates figures before compiling the book and independent chapters. `tools/build_figure_comparator.py` refreshes current numbers and captions and handles an empty pending list with a completion page. `tools/document_figure_project.py` updates source guides, index and decision ledger. `tools/audit_figure_project.py` verifies exact reading inputs, historical and data hashes, publication fonts, labels, dates and links.

The initial preparation/migration scripts preserve project history and must not be rerun over the current acceptance metadata. Acceptance is independent of the rendering engine; the manifest remains capable of tracking future revisions.

The final mesh correction uses quasiuniform triangles in smooth regions and thin contour-following strips, with explicit connectivity at the jump. Measured angles and aspect ratios verify tangential alignment. Both sides of the synthetic edge have smooth nonaffine intensities. The enhanced gradient display applies a common signed contrast curve to red/blue components without changing the numerical gradient or its magnitude.

See [final-integration.md](final-integration.md), `author-decisions.json`, `verification.json` and the per-figure guides for the correction record. Local build and visual-review evidence is retained under `build/final-figure-integration/`.
