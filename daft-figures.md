---
layout: default
body_class: figure-gallery-page
permalink: /daft/figures/
title: "DAFT figure gallery"
description: "All figures from The Discrete Algebra of the Fourier Transform, with searchable previews, vector PDFs, Python notebooks, and TikZ sources."
artwork: natures-mortes-2
figure_gallery: true
math: true
---
{% include art-banner.html %}
<main id="main" class="page-content wrap">
  <header class="page-heading gallery-heading">
    <p class="eyebrow"><a href="{{ '/daft/' | relative_url }}">← The DAFT book</a></p>
    <h1>A visual companion to <em>Fourier algebra</em></h1>
    <p class="page-lead">Explore every figure from <cite>The Discrete Algebra of the Fourier Transform</cite>. Find an idea, open it in the book, or use the code to make it your own.</p>
    <div class="gallery-intro-links"><span><strong>{{ site.data.daft_figures.figures.size }} figures</strong> · Python &amp; TikZ</span><a href="{{ site.data.daft_figures.archive | relative_url }}" download>Download all figure code and data <span aria-hidden="true">↓</span></a></div>
  </header>
  <div class="figure-tools" hidden>
    <div class="field"><label for="figure-search">Search the figures</label><input id="figure-search" type="search" placeholder="Try 3.2, Haar, interpolation, or heat" autocomplete="off" aria-controls="figure-grid"></div>
    <div class="field"><label for="figure-chapter">Chapter</label><select id="figure-chapter" aria-controls="figure-grid"><option value="">All chapters</option>{% for chapter in site.data.daft_figures.chapters %}<option value="{{ chapter.number }}">{{ chapter.number }} · {{ chapter.title }}</option>{% endfor %}</select></div>
    <div class="field"><label for="figure-type">Code</label><select id="figure-type" aria-controls="figure-grid"><option value="">Python &amp; TikZ</option><option value="Python">Python + notebook</option><option value="TikZ">TikZ</option></select></div>
    <div class="figure-filter-status"><p id="figure-count" role="status" aria-live="polite" aria-atomic="true">{{ site.data.daft_figures.figures.size }} figures</p><button type="button" id="reset-figures" hidden>Clear filters</button></div>
  </div>
  <noscript><p class="gallery-help">All figures are listed below in book order. Use your browser’s Find command to search this page.</p></noscript>
  <p class="gallery-help">Select a preview to enlarge it. Python notebooks open in Colab or download for Jupyter; TikZ diagrams include their LaTeX source.</p>
  <div class="figure-grid" id="figure-grid">
    {% for figure in site.data.daft_figures.figures %}
    <article class="figure-card" id="{{ figure.id }}" data-figure data-number="{{ figure.number }}" data-chapter="{{ figure.chapter }}" data-type="{{ figure.type }}" data-search="Figure {{ figure.number }} {{ figure.title | escape }} {{ figure.caption | escape }} {{ figure.chapter_title | escape }} {{ figure.type }}">
      <a class="figure-preview" href="{{ figure.pdf | relative_url }}" data-enlarge data-preview="{{ figure.preview | relative_url }}" aria-label="Enlarge Figure {{ figure.number }}: {{ figure.title | escape }}">
        <img src="{{ figure.thumbnail | relative_url }}" alt="{{ figure.title | escape }}" width="{{ figure.width }}" height="{{ figure.height }}" loading="lazy" decoding="async">
        <span class="figure-enlarge" aria-hidden="true">Enlarge ↗</span>
      </a>
      <div class="figure-card-copy">
        <div class="figure-meta"><a href="#{{ figure.id }}" class="figure-permalink">Figure {{ figure.number }}</a><span>{{ figure.type }}</span></div>
        <h2>{{ figure.title }}</h2>
        <p class="figure-caption">{{ figure.caption_html }}</p>
        <div class="figure-links">
          <a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page={{ figure.pdf_page }}">In the book <span class="figure-page">p. {{ figure.book_page }}</span></a>
          <a href="{{ figure.pdf | relative_url }}" download>PDF <span aria-hidden="true">↓</span></a>
          <a href="{{ figure.code }}">Python source <span aria-hidden="true">↗</span></a>
          {% if figure.tikz %}<a href="{{ figure.tikz }}">TikZ source <span aria-hidden="true">↗</span></a>{% endif %}
          {% if figure.notebook %}<a class="notebook-link" href="{{ figure.colab }}">Run notebook <span aria-hidden="true">↗</span></a><a href="{{ figure.notebook | relative_url }}" download>Notebook <span aria-hidden="true">↓</span></a>{% endif %}
        </div>
      </div>
    </article>
    {% endfor %}
  </div>
  <div class="empty-results" id="no-figures" hidden><p>No figures match these filters. Try a shorter term or choose another chapter.</p><button type="button" id="reset-empty-figures">Show all figures</button></div>
  <p class="gallery-footnote">Figures are numbered as in the book, including those in the exercise solutions. Chapter 8 has no numbered figures. <a href="{{ '/daft/' | relative_url }}">Return to the book</a>.</p>
</main>
<dialog id="figure-viewer" aria-labelledby="viewer-title">
  <div class="viewer-heading"><h2 id="viewer-title"></h2><button class="viewer-close" type="button" aria-label="Close enlarged figure">Close <span aria-hidden="true">×</span></button></div>
  <div class="viewer-image"><img id="viewer-image" alt=""><p id="viewer-image-error" hidden>The preview could not load. You can still open the figure PDF below.</p></div>
  <div class="viewer-footer"><div id="viewer-links" class="figure-links"></div><div class="viewer-navigation"><button id="previous-figure" type="button" aria-label="Previous figure">←</button><span id="viewer-position" aria-live="polite"></span><button id="next-figure" type="button" aria-label="Next figure">→</button></div></div>
</dialog>
