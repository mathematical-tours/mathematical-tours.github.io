# Interactive book

[Read the interactive edition online](https://mathematical-tours.github.io/book/interactive/).

The interactive reading edition is in `myst/`, using MyST and the same book theme as the author's [OT4ML reference](https://github.com/gpeyre/ot4ml/tree/main/myst). The canonical LaTeX and the published PDF editions remain intact.

- [Build, reading, and editing instructions](../../myst/README.md).
- [Figure-by-figure controls and exact publication captions](figures.md).
- [Validation summary](validation.json).

All 19 chapters, the preface, and the bibliography are converted from their existing sources. Equation, theorem, section, and figure numbering follow the PDF, including named mathematical tags. Browser-specific grouping of expanded accent macros preserves the intended mathematics. References to a merged figure resolve to its current figure; a label containing an apostrophe has a stable equivalent web anchor.

The 145 accepted numerical figures have interactive views across 30 mathematical families and 3,911 computed parameter combinations. The 81 accepted TikZ displays and 10 chapter illustrations explicitly retained by the author remain static; the original cover is preserved. Each numerical experiment retains the approved PDF figure as its comparison view.

Each experiment now supplies an exploration prompt, suggested comparisons, parameter explanations, and a short account of what to observe. Readers can inspect curves with the keyboard, keep numerical tables open while changing parameters, and open a shareable figure link that preserves their settings. Figure controls and plots follow the reading theme. Image changes are committed together only after the selected assets are available; failures show the approved illustration with a retry action.

Sliders select explicitly computed parameter values. The underlying routines are shared with the scientific figure generators, with fixed random samples and stated boundary conventions. This makes the complete book distributable as a static site without a running notebook service. The experiment notes distinguish schematic localization from actual curvelet computations, Fourier sampling from Radon projections, algebraic variance-stabilization inverses from unbiased ones, and finite Monte Carlo estimates from theorem statements.

Validation includes the complete MyST strict build; all local links and mathematical rendering; complete parameter grids and suggested comparisons that resolve to computed states; all figure panels in a browser; representative desktop and mobile layouts; and numerical invariants for transforms, approximation budgets, sensing adjoints, inpainting constraints, tangent-aligned meshes, Huffman coding, dual certificates, and Lloyd iterations. Additional checks cover histogram mass, the normalized logistic upper bound, curvelet translation directions, and unaliased birdsong detail. Reading tests block external requests and exercise search, section navigation, keyboard curve inspection, persistent numerical tables, saved links, comparison tabs, theme changes, rapid image updates, failed-request recovery, and responsive frame heights. Mathematical fonts and styles are included in the static output. Detailed local logs remain under `build/interactive-book/`.

The export shares one copy of each PDF and experiment asset tree. The unchanged PDF downloads retain the publication folder structure, and all 151 links from standalone chapters to named destinations in the full book are checked.

The GitHub Pages edition lives at `/book/interactive/`, linked from the site's existing book page. All 3,627 export files were checked byte for byte after a local Jekyll 3.9.5 build in safe mode. The browser audit uses that Jekyll output at the public path. Canonical page links, sharing images, the sitemap, and the PDF download menus use the publication address. The generated site and its checksum manifest are checked in under the parent repository's `book/interactive/` directory.
