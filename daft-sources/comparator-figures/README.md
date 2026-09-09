# Figure audit

Open [comparison.pdf](comparison.pdf): one original/revised pair per page,
preceded by an introductory page. The revised panels are the exact PDFs used
in the new book.

Each pair gives the **figure number, full caption, and printed book page** from
the compiled edition, together with the original figure number. Refer to that
number when requesting a change. `index.json` also identifies its stable label,
generated PDF and figure source.

There are 68 comparisons. Two original figures shared a single LaTeX figure
environment; they are listed separately here. The cover artwork is retained
and is not one of the numbered figures.

From the repository root, run `make comparator` after rebuilding the book.
The previous panels are reconstructed from the verified files under `legacy/`.
