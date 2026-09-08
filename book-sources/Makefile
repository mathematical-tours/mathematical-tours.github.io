.DEFAULT_GOAL := all
.PHONY: all book chapters figures figure-audit check clean help

PYTHON ?= python3
FIGURE_PYTHON ?= build/figure-runtime/bin/python
export FIGURE_PYTHON

all:
	$(PYTHON) tools/build_book.py
	$(MAKE) check
	$(PYTHON) tools/document_figure_project.py

book:
	$(PYTHON) tools/build_book.py --book-only

chapters: all

figures:
	$(FIGURE_PYTHON) tools/regenerate_figures.py

figure-audit check:
	$(FIGURE_PYTHON) tools/audit_figure_project.py > build/figure-regeneration/audit-console.txt
	@echo "Publication audit passed. Report: build/figure-regeneration/project-audit.json"

clean:
	$(PYTHON) tools/clean_build.py

help:
	@echo "make          Build and audit the book, all chapters, and figure review PDF"
	@echo "make book     Build the book and figure review PDF"
	@echo "make figures  Regenerate figure assets"
	@echo "make check    Audit the published book, chapters, and figures"
	@echo "make clean    Remove build caches; keep published PDFs and Python runtime"
	@echo "make web      Build and check the complete interactive MyST edition"
	@echo "make web-serve  Open a local static server for the interactive book"
	@echo "make web-publish  Prepare the validated edition for GitHub Pages"

# Interactive edition. Install its pinned Node dependencies once in myst/.
.PHONY: web web-sync web-data web-check web-serve web-browser-check web-document web-publish
web-sync:
	$(PYTHON) tools/interactive_book/convert.py

web-data: web-sync
	$(FIGURE_PYTHON) tools/interactive_book/build_data.py

web: web-data
	cd myst && npm run build -- --strict
	$(MAKE) web-check

web-check:
	$(PYTHON) tools/interactive_book/test_conversion.py
	$(FIGURE_PYTHON) tools/interactive_book/test_numerics.py
	cd myst && npm test
	$(PYTHON) tools/interactive_book/validate.py
	$(FIGURE_PYTHON) tools/interactive_book/validate_pdf_links.py

web-serve:
	cd myst && npm start

web-browser-check:
	cd myst && npm run test:browser

web-document:
	$(PYTHON) tools/interactive_book/document.py

# Copy a validated build to the existing GitHub Pages publication directory.
web-publish:
	$(PYTHON) tools/interactive_book/publish.py
