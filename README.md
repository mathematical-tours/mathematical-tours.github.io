# Mathematical Tours

[mathematical-tours.github.io](https://mathematical-tours.github.io/) presents three books by Gabriel Peyré:

- **Mathematical Foundations of Data Sciences:** PDF, separate chapters, and an interactive edition.
- **The Discrete Algebra of the Fourier Transform (DAFT):** an English adaptation of the French Ellipses book, with Python examples.
- **An Introduction to Data Sciences:** complete English and French books, and separate chapters.

The site also links to Numerical Tours and Mathematical Nexus and hosts the algorithms directory, course slides, and the ten-lecture *Mathematics of Machine Learning* course.

## Run the website locally

With Ruby and Bundler installed:

```sh
bundle install
python3 tools/prepare-site.py --destination .pages-source --report .pages-inputs.json
bundle exec jekyll serve --source .pages-source
```

Open `http://127.0.0.1:4000/`. The versioned `.github/workflows/pages.yml` workflow publishes changes to `master`. It stages only the pages and downloads declared in `tools/public-site.json`. The Gemfile follows GitHub Pages’ dependency bundle so local and published builds use the same versions. The templates also work with Jekyll 4.

For a static build and a check of the public pages and their local links:

```sh
bundle exec jekyll build --source .pages-source --destination _site
python3 tools/check-site.py _site
python3 tools/check-publication.py _site --source-map .pages-inputs.json
```

## Where to edit

| Content | Location |
| --- | --- |
| Homepage and book cards | `index.md`, `_data/books.json` |
| Book presentation pages | `book.md`, `daft.md`, `book-basics.md` |
| Algorithms and slide decks | `_data/algorithms.json`, `_data/slides.json` |
| Ten lectures and their materials | `_data/lectures.json` |
| Shared navigation, footer, metadata, and banners | `_includes/`, `_layouts/` |
| Website styles and navigation/search | `css/site.css`, `js/site.js` |
| Book typography and chapter catalogue | `css/book.css`, `js/book.js` |

The introduction page is generated from `book-basics-sources/website/`; update those templates when changing its layout or content. Book source folders have their own build instructions. The interactive edition is built from `book-sources/myst/` into `book/interactive/`; do not edit its generated HTML directly.

Published URLs remain stable, including `/maths-ia-course/`. The former Tweets pages forward readers to Mathematical Nexus. Development artifacts and archived sources are excluded from the website build; sources remain accessible on GitHub.

## Artwork and previews

The website banners and book covers use details from paintings by
[Louis Peyré (1923–2012)](https://www.louispeyre.com/). The original supplied
JPEGs live in `img/louis-peyre/`; `_data/artworks.json` assigns descriptive
alternative text and `_includes/art-banner.html` supplies the visible credit.
Choose a page's banner with its `artwork` front-matter key. The homepage's
“The paintings behind the pages” section explains their provenance.

Each book keeps a local copy of its cover artwork and a credit/checksum record,
so its PDF remains independently buildable. After changing a cover, rebuild the
book, its chapters and downloadable copies, then refresh the website thumbnail.
See `docs/artwork-publication.md` for the current cover assignments.

Resource previews show actual teaching material. Each entry in `_data/slides.json`
records its selected slide number, local JPEG, alternative text and original
Speaker Deck image URL. Preview links use `?slide=N` to open that exact slide.
`docs/visual-assets.json` records the preview sources and checksums.

## Deployment checks

The publishing workflow checks the site’s links, verifies every download and
interactive asset against its source checksum, and enforces a 600 MB / 5,000-file
budget before upload. Book sources, old datasets, caches, and unlinked auxiliary
materials stay in the repository without being copied to Pages. Existing chapter
URLs and compatibility downloads remain available.

To add a downloadable file, include it in `tools/public-site.json`; an undeclared
local link makes the build fail before publication. Preparation requires a fresh
destination directory, so it cannot overwrite unrelated local work.

The deployment step uses GitHub’s Pages REST API with the workflow’s temporary
OIDC identity and token. It waits for confirmed success, cancels on timeout or
interruption, and never logs credentials. This avoids the deprecated artifact SDK
in the current `deploy-pages` action. Official checkout and upload actions are
pinned to releases that support Node 24. See [the deployment repair record](docs/deployment-repair.md).
