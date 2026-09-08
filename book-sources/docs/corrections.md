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

## 8 September 2026: interactive content polishing

- Added mathematical exploration prompts, two or three suggested comparisons, and concise **What to look for** explanations to all 145 interactive figures. Checked each suggestion against its computed parameter grid. Parameter explanations distinguish wavelet dilation from decomposition depth, real coefficient budgets from encoded sizes, inverse regularization from penalty strength, and display gain from a change to the gradient operator.
- Refined figure typography, control grouping, discrete slider labels, legends, numerical tables, and responsive layouts. Kept boxed scientific axes, adapted plot labels and tick density to available width, and added a dark palette that follows the book's reading theme. All plotted series remain available in the legend, with long legends expandable.
- Added keyboard inspection of plotted samples: left/right select points, up/down select series, and Home/End reach endpoints. A marker and text readout identify the selection. Comparison tabs now support the same keyboard navigation conventions. Expanded numerical tables stay open after parameter changes.
- Added **Open figure** links that preserve the selected parameters and comparison tab, standalone figure headings with the publication number, and links back to the corresponding chapter figure. Image updates wait for the selected assets and ignore superseded slider requests. Failed data requests reveal the approved book illustration with a working retry action. The birdsong audio player is preserved across parameter changes.
- **Figure 1.1:** replaced direct waveform decimation with a peak-preserving overview, marked the selected interval in red, and displayed every recorded sample in the short-time zoom. The detail now ranges from approximately 80 ms to 5 ms, avoiding spurious oscillations from unfiltered subsampling.
- **Figure 5.24:** aligned the curvelet translation lattice with the selected frame atom: spacing is `1/sqrt(r)` tangentially and `1/r` normally. The spatial image uses the same upward-positive coordinate convention as the grid, and the actual discrete frequency orientation is reported.
- **Figures 12.2 and 12.3:** retained finite-sample eigenvalues outside the limiting Marchenko–Pastur support. Histograms now include every eigenvalue, normalize to unit mass, and draw bars at their true bin width. Added the number of eigenvalues outside the asymptotic support and clarified its finite-sample interpretation.
- **Figures 13.8 and 13.9:** included the energy of the final Lloyd iterate, so the current-iteration marker, centroids, and reported energy agree at every displayed step, including iteration 40.
- **Figure 13.16:** added the `logistic / log(2)` curve alongside the unscaled likelihood loss. The normalized curve touches the binary loss at zero margin and upper-bounds it everywhere. Included zero exactly in the plotting grid, exposed both normalizations in the current-value table, and stated the zero-margin convention.

Validation: 4 conversion tests, 15 numerical invariant tests, and 6 browser-module tests pass. The strict web build, all 236 figure displays, all 19 chapter pages, all 3,911 parameter states, and all 151 chapter-to-book PDF links pass their audits. Additional browser checks cover keyboard inspection, persistent tables, saved links, phone layouts, theme propagation, rapid parameter changes, and failed-request recovery. Approved book figures, canonical LaTeX, and published PDFs remain unchanged.

## 8 September 2026: book presentation page

- Redesigned the parent repository's `book.md` as a responsive HTML presentation page, with a dedicated Jekyll layout, local styles, and a cover preview rendered from the current full-book PDF.
- Added prominent links to the interactive edition and complete PDF, the author and affiliation, and the current book's page count and download size.
- Organized all 19 chapters in publication order into three topic groups. Each entry includes its canonical title, a short description, its PDF page count, and separate online/PDF links. Corrected outdated display titles, notably **Fourier and Convolution**, **Basics of Machine Learning**, **Shallow Learning**, and **Nonsmooth Convex Optimization**, while preserving the working PDF filenames.
- Added accessible chapter search and topic filtering. The complete chapter catalogue remains readable when JavaScript is disabled.
- Retained the Numerical Tours and the three optimization course-note PDFs. The Optimal Transport section now points to **OT4ML: Optimal Transport for Machine Learners**, at `https://www.gpeyre.com/ot4ml/`, as requested.
- Kept the established `/book/` page address and used a standalone layout so other website pages retain their existing presentation.

Validation: the full-book PDF, interactive homepage, all 19 chapter PDFs, all 19 web chapters, and the three course-note PDFs returned successful live responses. The generated page matches the canonical chapter order, titles, and PDF sizes. Browser checks passed at desktop, tablet, and 320-pixel mobile widths, including search, topic filters, keyboard operation, reset, chapter navigation, and the JavaScript-free catalogue. The page requires no external font or script downloads.
