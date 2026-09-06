# Completed figure review

All 145 reviewed reconstructions are integrated into the main book and its 19 independent chapters. There are **no pending comparisons**. [figure-comparisons.pdf](figure-comparisons.pdf) is now a one-page completion record.

The [previous 43-entry audit](archive/2026-09-06-reviewed-figures/figure-comparisons.pdf), its matching book snapshot and its numbering index are archived together. Its links still target that historical book. `author-decisions.json` records current and previous numbers for every inventory entry; `figure-project.json` records the accepted assets, source folders, exact captions and numerical provenance.

The 81 earlier accepted TikZ displays are retained. The original cover, Figures 11.3 and 11.7, and Figures 15.8–15.15 remain as requested. Old Figure 6.6 is removed; old Figures 7.9 and 7.10 are represented by the accepted combined Figure 7.9. The additional concrete Haar matrix is Figure 4.14.

## Organization

- Assets: `../../figures/chapter-<title>/<figure-name>/`.
- Sources: `../../figures-code/chapter-<title>/<figure-name>/`.
- Preserved historical inputs and composite: `original/` and `original.pdf`.
- Accepted reconstruction: `proposed.pdf`; the filename is retained for stable references.
- Numerical data, parameters and checks: `provenance.json`.

The default photograph is the supplied flower, with mandrill and a public-domain Felix frame used for the explicitly requested examples. No new reconstruction uses Lena. Historical panels remain in the archives.

## Build and verification

Run `make` from `book-sources/` to regenerate changed figures and rebuild the book, chapter PDFs and review status. Run `make figure-audit` to check the exact accepted reading inputs, original/data hashes, embedded fonts, labels, dates and PDF destinations. The latest build audit is `build/figure-regeneration/project-audit.json` (relative to `book-sources/`); `verification.json` records the latest publication verification.

The review generator still supports future pending entries: mark them explicitly in the manifest, then use the normal build or an unpublished `--preview` build. It extracts their exact current figure numbers and captions from the compiled book. Completed entries remain absent from the comparison index.

See [final-integration.md](final-integration.md) for this pass, [PROJECT.md](PROJECT.md) for the complete project policy, and the per-figure guides for mathematical details.
