---
layout: default
body_class: figure-gallery-page
permalink: /book-basics/figures/
title: "Introduction to Data Sciences · Figure gallery"
description: "All 42 illustrations from An Introduction to Data Sciences, in English and French, with searchable previews and PDF or JPG downloads."
artwork: paysages-2
gallery_book: basics
figure_gallery: true
math: true
---
{% include art-banner.html %}
<main id="main" class="page-content wrap">
  <header class="page-heading gallery-heading">
    <p class="eyebrow"><a href="{{ '/book-basics/' | relative_url }}">← An Introduction to Data Sciences</a></p>
    <h1>The ideas <em>in pictures</em></h1>
    <p class="page-lead">Explore every figure from <cite>An Introduction to Data Sciences</cite>. Browse the English or French illustrations, find them in the book, and download them for a closer look.</p>
    <div class="gallery-intro-links"><span><strong>42 figures</strong> · 5 chapters · English &amp; French</span><span>PDF &amp; JPG downloads</span></div>
  </header>
  <div class="figure-tools" hidden>
    <div class="field"><label for="figure-search">Search the figures</label><input id="figure-search" type="search" placeholder="Try entropy, noise, transport, or couleur" autocomplete="off" aria-controls="figure-grid"></div>
    <div class="field"><label for="figure-chapter">Chapter</label><select id="figure-chapter" aria-controls="figure-grid"><option value="">All chapters</option>{% for chapter in site.data.basics_figures.chapters %}<option value="{{ chapter.number }}">{{ chapter.number }} · {{ chapter.title }}</option>{% endfor %}</select></div>
    <div class="field"><label for="figure-type">Edition / Édition</label><select id="figure-type" data-query="edition" data-default="en" aria-controls="figure-grid"><option value="en">English</option><option value="fr">Français</option><option value="">Both editions / Les deux</option></select></div>
    <div class="figure-filter-status"><p id="figure-count" role="status" aria-live="polite" aria-atomic="true">42 figures</p><button type="button" id="reset-figures" hidden>Clear filters</button></div>
  </div>
  <noscript><p class="gallery-help">Both editions are listed below. Use your browser’s Find command to search in English or French.</p></noscript>
  <p class="gallery-help">Select a preview to enlarge it. Search works in both languages. / <span lang="fr">Cliquez sur une figure pour l’agrandir ; la recherche accepte les deux langues.</span></p>
  <div class="figure-grid" id="figure-grid">
    {% for figure in site.data.basics_figures.figures %}
    <article class="figure-card" id="{{ figure.id }}" lang="{{ figure.lang }}" data-figure data-number="{{ figure.number }}" data-chapter="{{ figure.chapter }}" data-type="{{ figure.lang }}" data-search="Figure {{ figure.number }} {{ figure.search | escape }}">
      <a class="figure-preview" href="{{ figure.pdf | relative_url }}" data-enlarge data-preview="{{ figure.jpg | relative_url }}" aria-label="{% if figure.lang == 'fr' %}Agrandir{% else %}Enlarge{% endif %} Figure {{ figure.number }}: {{ figure.title | escape }}">
        <img src="{{ figure.thumbnail | relative_url }}" alt="{{ figure.title | escape }}" width="{{ figure.width }}" height="{{ figure.height }}" loading="lazy" decoding="async">
        <span class="figure-enlarge" aria-hidden="true">{% if figure.lang == 'fr' %}Agrandir{% else %}Enlarge{% endif %} ↗</span>
      </a>
      <div class="figure-card-copy">
        <div class="figure-meta"><a href="#{{ figure.id }}" class="figure-permalink">Figure {{ figure.number }}</a><span>{% if figure.lang == 'fr' %}Français{% else %}English{% endif %}</span></div>
        <h2>{{ figure.title }}</h2>
        <p class="figure-caption">{{ figure.caption_html }}</p>
        <div class="figure-links">
          <a href="{{ figure.book | relative_url }}#page={{ figure.pdf_page }}" hreflang="{{ figure.lang }}">{% if figure.lang == 'fr' %}Dans le livre{% else %}In the book{% endif %} <span class="figure-page">p. {{ figure.book_page }}</span></a>
          <a href="{{ figure.pdf | relative_url }}" download="intro-data-sciences-figure-{{ figure.number }}-{{ figure.lang }}.pdf">PDF <span aria-hidden="true">↓</span></a>
          <a href="{{ figure.jpg | relative_url }}" download="intro-data-sciences-figure-{{ figure.number }}-{{ figure.lang }}.jpg">JPG <span aria-hidden="true">↓</span></a>
        </div>
      </div>
    </article>
    {% endfor %}
  </div>
  <div class="empty-results" id="no-figures" hidden><p>No figures match these filters. Try a shorter term or choose another chapter.</p><button type="button" id="reset-empty-figures">Show all figures</button></div>
  <p class="gallery-footnote">Each download contains the complete illustration from the selected edition, including panel labels. Captions retain the book’s figure and reference numbers. Flower photograph: Maïtine Bergounioux; other attributions appear in the captions and the book. <a href="{{ '/book-basics/' | relative_url }}">Return to the book</a>.</p>
</main>
<dialog id="figure-viewer" aria-labelledby="viewer-title">
  <div class="viewer-heading"><h2 id="viewer-title"></h2><button class="viewer-close" type="button" aria-label="Close enlarged figure">Close <span aria-hidden="true">×</span></button></div>
  <div class="viewer-image"><img id="viewer-image" alt=""><p id="viewer-image-error" hidden>The preview could not load. You can still open the figure PDF below.</p></div>
  <div class="viewer-footer"><div id="viewer-links" class="figure-links"></div><div class="viewer-navigation"><button id="previous-figure" type="button" aria-label="Previous figure">←</button><span id="viewer-position" aria-live="polite"></span><button id="next-figure" type="button" aria-label="Next figure">→</button></div></div>
</dialog>
