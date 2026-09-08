# Interactive edition

[Read the published book](https://mathematical-tours.github.io/book/interactive/).

The interactive edition follows the [OT4ML MyST book](https://www.gpeyre.com/ot4ml/myst/_build/html/) and uses the same MyST `book-theme`. It contains the preface, all 19 active chapters, the bibliography, and all 236 chapter figure displays. The 145 accepted numerical figures have interactive experiments; the 81 previously accepted TikZ displays and 10 expressly retained originals keep their approved appearance. The original cover is also retained.

The LaTeX in `../chapters-tex/` remains the source of truth. The web build does not modify the book or standalone PDFs. Their exact copies are available from the chapter pages. The download tree retains `FundationsDataScience.pdf` beside `chapters-pdf/`, preserving the chapter PDFs' relative links back to the book.

## Build and read

Requirements: Node.js 20 or later, Pandoc 3 (tested with 3.1.11), and the scientific Python environment described in the repository README. Node dependencies are pinned in `package-lock.json`; Python dependencies remain in `../figures-code/requirements.txt`.

From the repository root:

```sh
cd myst
npm ci
cd ..
make web
make web-serve
```

Open [the local book](http://127.0.0.1:8790/). `make web` converts the chapters, computes the experiments, builds MyST in strict mode, and checks mathematical invariants, local links, complete figure coverage, and unchanged PDF downloads. Generated output is `myst/_build/html/`; generated Markdown, images, parameter data, and build products are ignored by Git.

A web rebuild does not require TeX Live when the published PDFs and the checked-in reference snapshot are current. After an editorial change, build the PDF editions first with `make`, then run `make web`. The converter compares equation counters with the compiled labels and refreshes `reference-metadata.json` and `references-source.tex`. Those small snapshots let a clean checkout reproduce the same numbering and bibliography.

For a different static hosting prefix:

```sh
BASE_URL=/your/book/path make web
```

Publish the contents of `_build/html/` at that prefix. Use an empty `BASE_URL` for a site root. The site needs an ordinary static HTTP server; there is no Python/Jupyter service, account, API key, or external data request at reading time. Opening HTML through `file://` is not supported because browsers restrict local data fetching.

Set `SITE_URL=https://your-domain.example` when hosting on a different domain. It defaults to `https://mathematical-tours.github.io` and supplies canonical links, social preview URLs, and the sitemap.

### Publish on the existing website

The public address is <https://mathematical-tours.github.io/book/interactive/>. The parent repository publishes GitHub Pages from `master` using Jekyll. A checked-in static export in `../book/interactive/` (relative to `book-sources/`) supplies this edition, while `book.md` remains the site's book landing page.

From `book-sources/`, prepare the edition with:

```sh
BASE_URL=/book/interactive make web
make web-serve
```

The preview server uses the build's hosting prefix automatically. In a second terminal, run:

```sh
BOOK_URL=http://127.0.0.1:8790/book/interactive make web-browser-check
python3 tools/interactive_book/document.py --base-url https://mathematical-tours.github.io/book/interactive
make web-publish
```

`web-publish` checks that the static and browser audits match this exact build and its public prefix. It copies only the generated edition, records every file's checksum in `publication.json`, and removes only obsolete files owned by the previous publication. Commit the source changes and the parent repository's `book/interactive/` directory together, then push `master`. GitHub Pages handles publication through its existing build.

The parent `_config.yml` includes the theme's `_assets` and `_shared` directories and `_index` route files, which Jekyll otherwise omits, and excludes the unbuilt `book-sources/myst/` source tree. The main site keeps its existing Jekyll configuration and pages.

The theme's download menu serves the corresponding published PDF. The build removes obsolete cached exports and does not publish generated Markdown downloads, which Jekyll would otherwise convert into separate HTML pages. Canonical LaTeX and the web conversion code remain available in this repository.

## Experiments

Every experiment begins with a short exploration prompt, suggested comparisons under **Try**, and a **What to look for** explanation. Parameter notes distinguish such choices as wavelet dilation versus decomposition depth, or display contrast versus a change to the underlying operator. **Reset** restores the default settings; **Explore / Book figure** switches to the approved PDF illustration.

Native range inputs and menus work with the keyboard. Focus a plot and use left/right to inspect its samples, up/down to select a series, and Home/End to reach the endpoints. A visible marker and text readout identify the selected value. **Numerical values** opens a table that stays open when parameters change. Comparison tabs support the arrow keys and Home/End.

**Open figure** opens a standalone view with the current parameters and comparison tab encoded in its URL, ready to bookmark or share. That view includes the publication figure number and a link back to its chapter. Figures follow the book's light/dark theme, and plot colors adapt for contrast. A failed experiment request exposes the approved book figure and a retry button; rapid slider changes display only the most recently selected result.

The experiments use **computed discrete parameter grids**. Each slider stop loads the result for that exact parameter combination; it does not interpolate pictures. Heavy reconstruction algorithms therefore work on a static site, including when the server is offline. Noise samples, datasets, initializations, and measurement patterns are kept fixed within comparisons. Notes beneath each experiment state boundary conventions, finite-sample assumptions, solver settings, and any simplification relative to the publication figure.

The numerical routines reuse the approved figure code in `../figures-code/`. Important conventions are retained: boxed plot axes, dyadic coefficient separators, real coefficient budgets for approximation, signed red/blue gradients, actual P1 interpolation on tangent-aligned meshes, and real JPEG/JPEG-2000 encoding. Flower is the default photographic input; the authorized mandrill, the public-domain film cartoon, the recorded bird, Iris, and handwritten digits are available where relevant. No Lena imagery is used.

See the [figure-by-figure audit guide](../docs/interactive-book/figures.md) for precise numbers, captions, and controls.

## Editing and maintenance

| File or directory | Responsibility |
| --- | --- |
| `../chapters-tex/` | Canonical prose, mathematics, references, and figure order. |
| `interactions.json` | Explicit mapping of figure IDs to experiment families and views. |
| `../tools/interactive_book/convert.py` | Balanced LaTeX preprocessing and Pandoc-to-MyST conversion. |
| `../tools/interactive_book/experiments_*.py` | Parameter choices and scientific experiment definitions. |
| `../tools/interactive_book/numerics.py` | Shared data, image, curve, and serialization helpers. |
| `../tools/interactive_book/build_data.py` | Parameter-grid computation, caching, fallbacks, and panel generation. |
| `live/book-live.mjs`, `live/plot.mjs` | Browser controls, accessible views, and boxed scientific plots. |
| `live/guidance.mjs` | Mathematical exploration prompts, parameter explanations, and suggested comparisons. |
| `book-theme.css`, `live/book-live.css` | Reader and experiment typography/layout. |
| `book-frame.js` | Responsive iframe sizing and section navigation after route changes. |
| `scripts/build.mjs` | Strict-capable MyST build and portable iframe paths. |
| `../tools/interactive_book/serve.py` | Local preview at the build's configured URL prefix. |
| `../tools/interactive_book/publish.py` | Audited static publication and checksum manifest for GitHub Pages. |

Use `make web-sync` to convert only the text, or `make web-data` to refresh text and experiment data. For one experiment, run:

```sh
build/figure-runtime/bin/python tools/interactive_book/build_data.py --id approximation--triangulation-image
```

The cache includes experiment code, shared figure routines, and source data. A full data build prunes only unreferenced generated image files. It never removes the approved figure assets or original data.

## Validation

```sh
make web-check
make web-serve                 # keep this running in one terminal
make web-browser-check         # run from another terminal
```

The browser test checks every figure's controls, reset and fallback, all chapter pages, mobile layouts, navigation, search, and representative screenshots. It also exercises keyboard plot inspection, persistent tables, saved parameter links, theme propagation, rapid image updates, and recovery after an intentionally failed request. It blocks requests outside the local book, testing that reading does not require an external service. If Playwright has no browser installed, run `cd myst && npx playwright install chromium`; alternatively set `BROWSER_EXECUTABLE` to an installed Chromium/Chrome executable. Set `BOOK_URL` to test another local base URL.

Reports and screenshots are generated under `../build/interactive-book/`. The publication summary is in [docs/interactive-book](../docs/interactive-book/README.md).

Run `make web-document` after successful checks to regenerate the numbered audit guide and publication summary. The summary requires matching content fingerprints from the static audit and browser report, so an incomplete or stale run cannot be reported as a validated edition.

The build bundles the KaTeX stylesheet and fonts matching the theme's math renderer. Only those assets and their license are copied from the pinned KaTeX package; its JavaScript is not shipped separately. The theme's unused external notebook-widget styles are removed from both HTML and client bundles. Reading, mathematical typesetting, search, and the figure controls therefore require no CDN connection.

The static exporter consolidates identical PDF downloads and removes unreferenced duplicate copies of the generated experiment and asset trees. It only prunes its generated HTML output, never the publication PDFs or source data.
