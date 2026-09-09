# Louis Peyré artwork and teaching previews

Publication update: 9 September 2026.

The page banners and full-book covers feature details from paintings by
[Louis Peyré (1923–2012)](https://www.louispeyre.com/), supplied by Gabriel Peyré
in `img/louis-peyre/`. The source JPEGs are preserved byte for byte. Filenames
identify the supplied images; they are not asserted to be catalogue titles.
The homepage explains this provenance beside the supplied self-portrait.

| Book | Supplied cover image | Attribution record |
| --- | --- | --- |
| Mathematical Foundations of Data Sciences | `paysages-paris.jpg` | `book-sources/docs/cover-artwork.json` |
| The Discrete Algebra of the Fourier Transform | `natures-mortes-2.jpg` | `daft-sources/docs/cover-artwork.json` |
| An Introduction to Data Sciences, English and French | `paysages-2.jpg` | `book-basics-sources/figures/chapter-book-cover/sunset-landscape/credit.json` |

Each cover has a visible artist credit linked to the artist's website. Full-book
cover images are copied into their respective source repositories so that a
downloaded book source tree does not depend on its parent website directory.
The archived books and former cover assets are preserved.

All four full books and 29 independent chapters were recompiled. The book and
figure-comparison checks verify embedded fonts, links, numbering and clean TeX
logs. The interactive MFDS edition uses the same Paris artwork and credits it on
the landing page and chapter pages. Its downloadable PDFs match the new builds.

The four companion-resource cards, the ten-lecture course card, and all 17 course
slide decks now have previews from the teaching material. The slide images are
selected interior slides, not repeated title slides; each preview opens its exact
slide using Speaker Deck's `?slide=N` URL. Descriptive alternatives, source URLs,
slide numbers and image checksums are recorded in `_data/slides.json` and
`docs/visual-assets.json`. These teaching previews are separate from the credited
painting banners. The preview selection avoids Lena imagery.
