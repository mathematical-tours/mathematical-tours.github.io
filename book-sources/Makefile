.PHONY: all book chapters figures figure-audit

FIGURE_PYTHON ?= build/figure-runtime/bin/python
export FIGURE_PYTHON

all:
	python3 scripts/build_book.py

book:
	python3 scripts/build_book.py --book-only

chapters: all

figures:
	$(FIGURE_PYTHON) scripts/regenerate_figures.py

figure-audit:
	$(FIGURE_PYTHON) scripts/audit_figure_project.py
