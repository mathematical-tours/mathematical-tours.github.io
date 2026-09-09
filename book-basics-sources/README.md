# An Introduction to Data Sciences

**Gabriel Peyré — CNRS & DMA, École Normale Supérieure, PSL**

English and French editions of an illustrated introduction to information theory,
image processing, sparse approximation and compressed sensing, optimal transport,
and neural networks. Both complete books contain the same five chapters.

- [Complete English book](IntroImaging.pdf)
- [Livre complet en français](IntroImaging-FR.pdf)
- [English chapters](chapters-pdf/en/) · [Chapitres en français](chapters-pdf/fr/)
- [Presentation page](https://mathematical-tours.github.io/book-basics/)
- [Searchable figure gallery](https://mathematical-tours.github.io/book-basics/figures/)
- [Major corrections](corrections.md)

Les deux éditions ont été relues et harmonisées. Les sources françaises sont en
UTF-8 normalisé ; les PDF utilisent des polices incorporées et permettent la copie
du texte accentué. Chaque chapitre possède une introduction, sa bibliographie et,
dans sa version autonome, le titre du livre, l'auteur, l'affiliation et la date.

## Organization

| Path | Contents |
| --- | --- |
| `IntroImaging.tex`, `IntroImaging-FR.tex` | Complete-book drivers; published PDFs use the same stems |
| `chapters-tex/en/`, `chapters-tex/fr/` | Chapter text, including the two neural-network sections |
| `chapters-tex/shared/` | Typography, metadata, macros and the common bibliography |
| `chapters-tex/standalone/` | Drivers for the ten independent chapter PDFs |
| `chapters-pdf/en/`, `chapters-pdf/fr/` | Published chapter downloads |
| `figures/` | Only the illustration assets used by the revised books |
| `gallery/` | Complete figure PDF/JPG downloads and lightweight previews for both editions |
| `data/`, `figures-code/` | Flower photograph and reproducible numerical illustrations |
| `tools/` | Builds, mathematical checks, PDF checks and page generation |
| `website/` | Bilingual Jekyll page template and its layout |
| `reports/` | Results of the latest verification and figure generation |
| `legacy/pre-revision-2026-09-08/` | Unmodified original sources, PDFs and associated material |
| `build/` | Local intermediate output, ignored by Git and Jekyll |

## Build the books

Install a current TeX Live or MacTeX distribution containing `latexmk`, pdfLaTeX,
BibTeX, `newpx`, `sourcesanspro`, `microtype`, `babel-french`, `titlesec` and TikZ.
The PDF build uses Python 3.10 or later and its standard library; no Python
packages are needed merely to compile the existing figures.

```sh
python3 tools/build.py
```

This produces both complete books and all ten chapter PDFs. TeX auxiliary files
are kept in separate `build/latex/` directories. The build fails if the final log
contains warnings, missing glyphs, undefined references, or overflowing boxes.
Use `--books` or `--chapters` to build a subset, and `--jobs 1` for a sequential build.
Run commands from this directory; the Python tools also accept other working directories.

## Reproduce and verify

Use a virtual environment for the numerical examples and PDF checks:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python figures-code/chapter-sampling/generate.py
.venv/bin/python tools/build.py
.venv/bin/python tools/check-mathematics.py
.venv/bin/python tools/check-pdfs.py
```

The mathematical checks cover the coding example, entropy calculations, image
transformations, transport assignment and its dual certificate, wavelet identities,
and neural-network gradients. The generator separately verifies its wavelet
approximation errors and the convergence residual of the inpainting computation.
The PDF audit checks all twelve files, embedded fonts, French text extraction,
chapter and figure numbering, and the placement of the edition date. Numerical
checks supplement the mathematical review; they are not substitutes for proofs.

The flower photograph comes from Maïtine Bergounioux's material in the original
book. Its source is retained in `data/flower.png`. The revised sampling examples
use this photograph throughout. The historical illustrations remain attributed
in their captions and the text; see [figure notes](figures/README.md).

## Presentation page

```sh
.venv/bin/python tools/build-website.py --output build/site
```

This prepares `book-basics.md`, `_layouts/book-basics.html` and
`img/book-basics-cover.webp`, using the actual page counts and file sizes of the
PDFs. The page reuses the parent website's `css/book.css`. To update the parent
Jekyll site, run the command with `--output ..`; this writes only those three
presentation files. It does not publish the site.

The complete-book URLs are unchanged. Older transport and neural-network download
URLs are maintained by the compatibility copies described in
`tools/legacy-downloads.json`. These copies contain the revised chapter PDFs and
are refreshed by the build.

Before changing the edition date, edit `chapters-tex/shared/book-metadata.tex`, then
rebuild and regenerate the page. The complete books show the date on the cover;
standalone chapters show it below their title.

## Figure gallery

The [figure gallery](https://mathematical-tours.github.io/book-basics/figures/)
contains all 42 numbered illustrations from each edition. Search accepts English
and French keywords, figure numbers, and chapter filters. Each entry offers a
larger preview, a link to the correct book page, and PDF and JPG downloads.

After changing a book or an illustration, rebuild the complete books and exports:

```sh
.venv/bin/python tools/build.py --books
.venv/bin/python tools/build-gallery.py
python3 ../tools/check-basics-gallery.py
```

The exporter uses pdfLaTeX (including the `preview` package) and the PDF/image
packages in `requirements.txt`. It reads captions, numbering and page destinations
from the current sources and compiled books. Complete figure bodies retain their
panels and labels; the PDFs preserve the original vector and raster content.
It also generates high-quality JPGs and smaller WebP previews. Visually identical
English/French illustrations share files; distinct labels keep separate exports.

By default, the tool writes `gallery/` and the parent site's
`_data/basics_figures.json`. Use `--site-root build/gallery-preview` to export into
a temporary site directory. The site build checks source and download hashes so
that changed illustrations, broken downloads or stale book links fail publication.
