---
layout: default
title: "Mathematical Tours — Books for Data Sciences"
description: "Three open books by Gabriel Peyré on the mathematics of data sciences, with Python examples, interactive figures, lectures, and computational resources."
body_class: home-page
---
{% include art-banner.html %}
<main id="main">
  <section class="home-intro wrap" aria-labelledby="home-title">
    <div><p class="eyebrow">Open books by Gabriel Peyré</p><h1 id="home-title">Mathematics for<br><em>data sciences.</em></h1></div>
    <div class="home-intro-copy"><p>Understand the ideas behind the methods. Three books connect mathematical theory with the study of signals, images, and learning.</p><a class="text-link" href="#books">Find your next book <span aria-hidden="true">↓</span></a></div>
  </section>
  <section class="bookshelf wrap" id="books" aria-labelledby="books-title">
    <div class="section-heading"><h2 id="books-title">The books</h2><p>Read a complete book or choose a chapter.<br>All editions are freely available.</p></div>
    <div class="book-grid">
    {% for book in site.data.books %}
      <article class="book-card book-card-{{ book.id }}">
        <a class="book-card-cover" href="{{ book.url | relative_url }}" tabindex="-1" aria-hidden="true"><img src="{{ '/img/' | append: book.image | relative_url }}" alt="" width="849" height="1200" loading="lazy"></a>
        <div class="book-card-copy"><p class="eyebrow">{{ book.tag }}</p><h3><a href="{{ book.url | relative_url }}">{{ book.title }}</a></h3><p class="book-summary">{{ book.summary }}</p><p class="book-level">{{ book.level }}</p><p class="book-formats">{{ book.formats }}</p><div class="book-card-actions"><a class="text-link" href="{{ book.url | relative_url }}">{{ book.action }} <span aria-hidden="true">→</span></a><a href="{{ book.pdf | relative_url }}" aria-label="{{ book.title }}: full PDF">PDF ↓</a></div></div>
      </article>
    {% endfor %}
    </div>
  </section>
  <section class="home-resources" id="resources" aria-labelledby="resources-title"><div class="wrap">
    <div class="section-heading"><div><p class="eyebrow">Beyond the books</p><h2 id="resources-title">Explore, compute, and teach</h2></div><p>Companion resources to see the ideas in action and follow them further.</p></div>
    <div class="companion-grid">
      <article class="companion-card"><span class="resource-index" aria-hidden="true">01 / COMPUTE</span><h3>Numerical Tours</h3><p>Learn through computational experiments: Python notebooks on images, optimization, geometry, and machine learning.</p><a class="text-link" href="https://www.numerical-tours.com/">Explore the notebooks <span aria-hidden="true">↗</span></a></article>
      <article class="companion-card"><span class="resource-index" aria-hidden="true">02 / DISCOVER</span><h3>Mathematical Nexus</h3><p>A searchable collection of interactive notebooks and visual explanations, bringing mathematical ideas into view.</p><a class="text-link" href="https://www.gpeyre.com/mathematical-nexus/">Visit the collection <span aria-hidden="true">↗</span></a></article>
      <article class="companion-card"><span class="resource-index" aria-hidden="true">03 / FIND A METHOD</span><h3>Algorithms</h3><p>Browse algorithms and mathematical models by subject, from coding and graphs to numerical analysis and learning.</p><a class="text-link" href="{{ '/algorithms/' | relative_url }}">Browse the directory <span aria-hidden="true">→</span></a></article>
      <article class="companion-card"><span class="resource-index" aria-hidden="true">04 / TEACH & REVIEW</span><h3>Course slides</h3><p>Two teaching collections covering signal and image processing, and the geometry of surfaces and meshes.</p><a class="text-link" href="{{ '/slides/' | relative_url }}">Find a slide deck <span aria-hidden="true">→</span></a></article>
    </div>
    <article class="course-feature"><div class="course-feature-number" aria-hidden="true">10</div><div><p class="eyebrow">A guided course · Ten two-hour lectures</p><h3>Mathematics of Machine Learning</h3><p>From gradients and sparse models to neural networks, generative models, and optimal transport. Follow the lectures with transcripts, notebooks, and slides.</p></div><a class="button primary" href="{{ '/maths-ia-course/' | relative_url }}">Explore the course <span aria-hidden="true">→</span></a></article>
  </div></section>
</main>
