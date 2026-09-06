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
