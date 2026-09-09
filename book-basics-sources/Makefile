PYTHON ?= python3
SITE_DIR ?= build/site

.PHONY: all pdf books chapters figures check website gallery clean
all: pdf
pdf:
	$(PYTHON) tools/build.py
books:
	$(PYTHON) tools/build.py --books
chapters:
	$(PYTHON) tools/build.py --chapters
figures:
	$(PYTHON) figures-code/chapter-sampling/generate.py
check: pdf
	$(PYTHON) tools/check-mathematics.py
	$(PYTHON) tools/check-pdfs.py
gallery: books
	$(PYTHON) tools/build-gallery.py
website:
	$(PYTHON) tools/build-website.py --output $(SITE_DIR)
clean:
	$(PYTHON) -c "import shutil; shutil.rmtree('build', ignore_errors=True)"
