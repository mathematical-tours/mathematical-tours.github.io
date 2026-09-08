# Kraft's inequality: detailed proof draft

- [kraft-inequality.pdf](kraft-inequality.pdf): four-page self-contained explanation and proof, with a worked example and a diagram.
- [kraft-inequality.tex](kraft-inequality.tex): editable standalone source.
- [book-insert.tex](book-insert.tex): a shorter replacement for the converse paragraph in the current book proof.

The current argument in `chapters-tex/shannon-coding.tex` is correct. Its compressed step is the alignment claim: for sorted lengths, every earlier block size is an integer multiple of the next size. The draft proves the exact correspondence between aligned leaf blocks and binary-tree vertices, constructs the codewords explicitly, and establishes the prefix property and the tree representation. It also handles the empty word and gives an alternative proof by counting available vertices.

The main book and its published PDFs are unchanged. This draft is intentionally separate from the active chapter includes.

Build from `book-sources/`:

```sh
mkdir -p build/draft-kraft
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build/draft-kraft draft/kraft-inequality/kraft-inequality.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build/draft-kraft draft/kraft-inequality/kraft-inequality.tex
cp build/draft-kraft/kraft-inequality.pdf draft/kraft-inequality/kraft-inequality.pdf
```

The proof uses the finite binary-alphabet setting and notation of Lemma `lem-kraft` in the book. No additional coding theorem is needed.

Validation: the standalone draft compiles without warnings or box diagnostics. The shorter fragment also compiles independently with the book preamble, and the draft pages have been visually checked.
