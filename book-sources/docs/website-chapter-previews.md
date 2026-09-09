# Illustrated chapter navigation

The Foundations presentation page now gives every chapter a representative
preview beside its description, online reading link and PDF download. Keyword
and topic filters continue to apply to the complete chapter entries. The layout
shares its responsive styles with the introductory book page.

Eighteen previews reproduce figures from the current book. Their image sources
and figure numbers are resolved from the published reader’s figure catalogue.
When a figure contains separate panel PDFs, the preview uses a complete panel.
Each thumbnail links to that figure’s actual anchor in the online chapter.

Chapter 16 has no numbered illustration. Its small SVG diagram represents the
model defined in `chapters-tex/perceptrons.tex`: two input coordinates, four
nonlinear hidden neurons, and a linear output combining their activations. It
links to that model’s equation in the reader. This illustration is for the
presentation page; the book text and PDFs are unchanged.

To regenerate the previews using the existing PDF/image environment:

```sh
python3 tools/build_chapter_previews.py
```

The command requires `pypdfium2` and `Pillow`. It writes
`../img/previews/foundations-*` and
`../_data/foundations_chapter_previews.json`. Use `--site-root build/site-preview`
for an isolated export. The metadata records the selected source, its checksum,
the descriptive alternative text, and the verified reader destination.

The previews reuse the book’s illustrations and require about 520 KiB in total.
The presentation page retains its Louis Peyré painting banner and credits.
