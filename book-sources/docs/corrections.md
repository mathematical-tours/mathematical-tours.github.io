# Editorial corrections

## 7 September 2026: prefix coding and Shannon's lower bound

Updated `chapters-tex/shannon-coding.tex` in the introduction to prefix codes and in the proof of Shannon's source coding theorem.

- Defined unique decodability and added the two-symbol code `0`, `01`, which is uniquely decodable even though `0` is a prefix of `01`. Explained the unique parsing of `00101` as `0 | 01 | 01`.
- Distinguished immediate prefix decoding from delayed decisions. A decoder that commits prematurely to `0` must backtrack when the next bit is `1`; this particular code can instead be decoded with one bit of lookahead. The text describes that overhead without claiming an unavoidable increase in asymptotic decoding complexity.
- Expanded the proof of the lower bound using `q_k = 2^{-ell_k}/Z` and `ell_k = -log_2(q_k) - log_2(Z)`. Displayed the substitution into `L(c) - H(p)` and the normalization of `p`.
- Kept the lower-bound derivation together with the explanation of its inequalities, so an intervening figure cannot separate them.
- Retained equation label `eq-proof-shannon` and its number (3.3), with the explicit chain `KL(p | q) - log_2(Z) >= KL(p | q) >= 0`. Explained each inequality and the equality conditions.

The separate [detailed Kraft proof draft](../draft/kraft-inequality/README.md) remains outside the main book.

Validation: rebuilt the 277-page book, all 19 standalone chapters, and the figure review PDF with zero diagnostics. The Shannon coding chapter now has nine pages. Equation (3.3) is preserved in both editions; its derivation and explanation remain together. Visually checked the edited pages and passed the publication audit, including all 767 shared-label checks.

## 8 September 2026: interactive reading edition

- Added a reproducible MyST edition of the preface, all 19 chapters, and bibliography, with the PDF's mathematical and figure numbering. The LaTeX remains canonical.
- Added interactive views for all 145 accepted numerical figures across 30 experiment families and 3,911 computed parameter combinations, with the approved book figure retained in a comparison tab. Existing TikZ figures, the author's retained originals, and the original cover are preserved.
- Reused the scientific figure routines for parameter sweeps, with fixed noise and measurement patterns, meaningful menus, and exact computed slider values. Preserved wavelet scale separators, boxed axes, real coefficient budgets, and the approved photographic data.
- Refined tangent-aligned P1 mesh experiments, including a boundary zoom and an explicitly encoded equal-budget JPEG-2000 comparison. Enhanced signed gradient displays and added pixel-gradient arrows.
- Preserved very small optimization errors with significant-digit export so logarithmic plots remain correct. Removed view-specific controls that had no visual effect.
- Adapted expanded TeX accent arguments and named equation tags for browser mathematics without changing the manuscript. Resolved the merged convexity figure's alias and the apostrophe-containing equation label in web links.
- Added keyboard-operable controls, numerical plot tables, responsive figure heights, chapter/section navigation, search, and exact PDF downloads.
- Bundled mathematical styles and fonts locally and removed unused external widget styles, so reading and search work with external requests blocked. Checked hosting beneath a subdirectory as well as at the site root.
- Polished light and dark reading themes, collapsed the initial mobile section outline, and gave the mobile search button an accessible name. Checked embedded resizing after parameter changes, opening numerical tables, and switching comparison views.
- Added readable scientific tick spacing, preserved valid axes for constant negative curves and very small positive values, and made plot readouts follow the nearest displayed curve. Signed gradient legends now explicitly distinguish negative, zero, and positive values.
- Consolidated duplicate PDF downloads and unused copied asset trees in the static export. Retained the original book/standalone-chapter folder relationship and verified all 151 remote PDF links against their actual named destinations.
- Prepared publication at `https://mathematical-tours.github.io/book/interactive/`, linked from the existing book page. Updated that page's chapter list and the auxiliary PDF paths affected by the earlier repository cleanup.
- Preserved the MyST theme assets and homepage routes through GitHub Pages' Jekyll build. Recreated generated export directories to remove obsolete cached files, and replaced the theme's generated-Markdown downloads with the corresponding publication PDFs.
- Added a preview server that respects the public URL prefix and a publication command that requires matching successful audits and records every published file's checksum.
- Set canonical page links, sharing images, and sitemap URLs to the public address. Verified all 3,604 export files byte for byte after Jekyll processing and checked the corresponding PDF in every chapter's download menu.

See the [interactive edition record](interactive-book/README.md) for the figure guide, build instructions, and validation results.
