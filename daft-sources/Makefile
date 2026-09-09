PYTHON ?= python3
LATEXMK ?= latexmk

.PHONY: all book listings figures comparator check documents python-download gallery clean
all: book gallery python-download

listings:
	$(PYTHON) tools/sync_listings.py

book: listings
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error DAFT-EN.tex

figures:
	$(PYTHON) tools/regenerate_figures.py

comparator: book
	$(PYTHON) tools/build_comparator.py

check:
	$(PYTHON) -W error tools/verify_programs.py

documents:
	$(PYTHON) tools/verify_documents.py

gallery: book
	$(PYTHON) tools/build_gallery.py

python-download:
	$(PYTHON) tools/package_python.py

clean:
	$(LATEXMK) -c DAFT-EN.tex
	$(PYTHON) -c "from pathlib import Path; [(Path('comparator-figures') / ('comparison.' + ext)).unlink(missing_ok=True) for ext in ('aux', 'log', 'out')]"
