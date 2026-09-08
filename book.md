---
layout: book
permalink: /book/
title: "Mathematical Foundations of Data Sciences"
description: "Read Gabriel Peyré’s book in PDF or as an interactive edition. Explore signal processing, inverse problems, optimization, and machine learning across 19 chapters."
book_pdf: /book-sources/FundationsDataScience.pdf
book_pages: 277
book_pdf_mb: 16.3
groups:
  - id: signals
    title: "Signal and image processing"
    range: "01–07"
  - id: inverse
    title: "Inverse problems and sparsity"
    range: "08–12"
  - id: learning
    title: "Learning and optimization"
    range: "13–19"
chapters:
  - number: 1
    slug: shannon
    title: "Shannon Sampling Theory"
    summary: "Sampling, aliasing, and the reconstruction of continuous signals."
    group: signals
    pages: 9
  - number: 2
    slug: fourier
    title: "Fourier and Convolution"
    summary: "Fourier representations, convolution, filtering, and the discrete transform."
    group: signals
    pages: 24
  - number: 3
    slug: shannon-coding
    title: "Shannon Coding Theory"
    summary: "Entropy, prefix codes, Huffman coding, and block coding."
    group: signals
    pages: 9
  - number: 4
    slug: wavelets
    title: "Wavelets"
    summary: "Multiresolution analysis, filter banks, and wavelet representations."
    group: signals
    pages: 25
  - number: 5
    slug: approximation
    title: "Linear and Nonlinear Approximation"
    summary: "Coefficient selection, approximation rates, and geometric image models."
    group: signals
    pages: 22
  - number: 6
    slug: compression
    title: "Compression"
    summary: "Quantization, source models, and image compression."
    group: signals
    pages: 7
  - number: 7
    slug: denoising
    title: "Denoising"
    summary: "Noise models, filtering, and wavelet thresholding."
    group: signals
    pages: 25
  - number: 8
    slug: variational-priors
    title: "Variational Priors and Regularization"
    summary: "Image gradients, diffusion, and total-variation regularization."
    group: inverse
    pages: 12
  - number: 9
    slug: inverse-problems
    title: "Inverse Problems"
    summary: "Observation operators, reconstruction, and ill-posed problems."
    group: inverse
    pages: 18
  - number: 10
    slug: sparse-regularization
    title: "Sparse Regularization"
    summary: "Sparse models, regularization paths, and iterative reconstruction."
    group: inverse
    pages: 11
  - number: 11
    slug: sparse-theory
    title: "Theory of Sparse Regularization"
    summary: "Dual certificates, recovery conditions, and stability."
    group: inverse
    pages: 15
  - number: 12
    slug: compressed-sensing
    title: "Compressed Sensing"
    summary: "Random measurements and reconstruction from incomplete data."
    group: inverse
    pages: 12
  - number: 13
    slug: machine-learning
    title: "Basics of Machine Learning"
    summary: "Principal components, clustering, regression, and classification."
    group: learning
    pages: 27
  - number: 14
    slug: optim-ml-smooth
    title: "Optimization & Machine Learning: Smooth Optimization"
    summary: "Gradient descent, smoothness, convexity, and convergence."
    group: learning
    pages: 18
  - number: 15
    slug: optim-ml-advanced
    title: "Optimization & Machine Learning: Advanced Topics"
    summary: "Stochastic, accelerated, and second-order optimization methods."
    group: learning
    pages: 17
  - number: 16
    slug: perceptrons
    title: "Shallow Learning"
    summary: "Perceptrons, shallow neural networks, and approximation."
    group: learning
    pages: 9
  - number: 17
    slug: deep-learning
    title: "Deep Learning"
    summary: "Deep networks, backpropagation, and learned representations."
    group: learning
    pages: 6
  - number: 18
    slug: convex-analysis
    title: "Convex Analysis"
    summary: "Convex sets and functions, subgradients, and duality."
    group: learning
    pages: 10
  - number: 19
    slug: optim-nonsmooth
    title: "Nonsmooth Convex Optimization"
    summary: "Proximity operators, splitting methods, and nonsmooth objectives."
    group: learning
    pages: 13
---
<main id="main">
  <section class="hero wrap" aria-labelledby="book-title">
    <div class="hero-copy">
      <p class="eyebrow">Signal processing · Optimization · Learning</p>
      <h1 id="book-title">Mathematical Foundations of <em>Data Sciences</em></h1>
      <p class="byline">Gabriel Peyré<span class="affiliation">CNRS &amp; DMA · École Normale Supérieure</span></p>
      <p class="hero-intro">A mathematical route from signals and images to inverse problems, optimization, and machine learning. Read the theory, explore the figures, and connect the ideas to computational methods.</p>
      <div class="hero-actions">
        <a class="button primary" href="{{ '/book/interactive/' | relative_url }}">Read the interactive book <span aria-hidden="true">↗</span></a>
        <a class="button" href="{{ page.book_pdf | relative_url }}" download="FundationsDataScience.pdf">Download the full PDF <span aria-hidden="true">↓</span></a>
      </div>
      <p class="format-note">Free to read · Full PDF: {{ page.book_pdf_mb }} MB</p>
    </div>
    <figure class="cover">
      <a href="{{ page.book_pdf | relative_url }}" aria-label="Open the complete book PDF"><img src="{{ '/img/foundations-cover.webp' | relative_url }}" alt="Cover of Mathematical Foundations of Data Sciences by Gabriel Peyré, with the original wave illustration" width="849" height="1200" fetchpriority="high"></a>
      <figcaption>The complete book<span>PDF edition · {{ page.book_pages }} pages</span></figcaption>
    </figure>
  </section>

  <div class="book-facts wrap" aria-label="The book at a glance">
    <div class="fact"><strong>19</strong><span>chapters, from foundations to applications</span></div>
    <div class="fact"><strong>145</strong><span>interactive figures to explore</span></div>
    <div class="fact"><strong>{{ page.book_pages }}</strong><span>pages in the complete PDF</span></div>
  </div>

  <section class="chapters wrap" id="chapters" aria-labelledby="chapters-title">
    <div class="section-heading">
      <div><p class="eyebrow">Contents</p><h2 id="chapters-title">Read chapter by chapter</h2></div>
      <p>Follow the book in order or go straight to a topic. Every chapter is available online and as a separate PDF.</p>
    </div>
    <div class="chapter-tools" hidden>
      <div class="field"><label for="chapter-search">Find a chapter</label><input type="search" id="chapter-search" placeholder="Try wavelets, sparsity, or learning" autocomplete="off" aria-controls="chapter-catalogue"></div>
      <div class="field"><label for="chapter-topic">Topic</label><select id="chapter-topic" aria-controls="chapter-catalogue"><option value="">All topics</option>{% for group in page.groups %}<option value="{{ group.id }}">{{ group.title }}</option>{% endfor %}</select></div>
      <p class="chapter-count" id="chapter-count" role="status" aria-live="polite">19 chapters</p>
    </div>
    <div id="chapter-catalogue">
      {% for group in page.groups %}
      <section class="chapter-group" aria-labelledby="heading-{{ group.id }}">
        <div class="group-heading"><h3 id="heading-{{ group.id }}">{{ group.title }}</h3><span>Chapters {{ group.range }}</span></div>
        {% assign group_chapters = page.chapters | where: 'group', group.id %}
        <ol class="chapter-list" start="{{ group_chapters.first.number }}">
          {% for chapter in group_chapters %}
          {% assign online_path = '/book/interactive/' | append: chapter.slug | append: '/' %}
          {% assign pdf_path = '/book-sources/chapters-pdf/' | append: chapter.slug | append: '.pdf' %}
          <li data-chapter="{{ chapter.slug }}" data-topic="{{ group.id }}" data-search="{{ chapter.number }} {{ chapter.title | escape }} {{ chapter.summary | escape }} {{ chapter.slug }}">
            <span class="chapter-number" aria-hidden="true">{% if chapter.number < 10 %}0{% endif %}{{ chapter.number }}</span>
            <div><h4 class="chapter-title"><a href="{{ online_path | relative_url }}">{{ chapter.title | escape }}</a></h4><p class="chapter-description">{{ chapter.summary }}</p></div>
            <div class="chapter-actions"><a href="{{ online_path | relative_url }}" aria-label="Read chapter {{ chapter.number }} online: {{ chapter.title | escape }}">Read online <span aria-hidden="true">↗</span></a><a class="chapter-pdf" href="{{ pdf_path | relative_url }}" aria-label="Chapter {{ chapter.number }} PDF: {{ chapter.title | escape }}, {{ chapter.pages }} pages">PDF <span class="pages">{{ chapter.pages }} pp.</span><span aria-hidden="true"> ↓</span></a></div>
          </li>
          {% endfor %}
        </ol>
      </section>
      {% endfor %}
    </div>
    <p class="no-results" id="no-chapters" hidden>No chapters match these filters. <button id="reset-chapters" type="button">Show all chapters</button></p>
  </section>

  <section class="about wrap" aria-labelledby="about-title">
    <h2 id="about-title">From mathematical ideas<br>to numerical methods</h2>
    <div><p>This book develops the mathematical foundations behind signal and image processing, inverse problems, and machine learning. Linear operators, nonlinear approximation, and convex optimization provide a common language for understanding the methods and turning them into algorithms.</p><p>It is written as a companion to the <a href="https://www.numerical-tours.com/">Numerical Tours of Data Sciences</a>. The interactive edition follows the same chapters and lets you change parameters, compare methods, and inspect the figures alongside the mathematics.</p></div>
  </section>

  <section class="resources" id="resources" aria-labelledby="resources-title">
    <div class="wrap">
      <p class="eyebrow">Companion material</p><h2 id="resources-title">Continue exploring</h2>
      <div class="resource-grid">
        <article class="resource-card"><span class="resource-type">Implementations</span><h3>Numerical Tours</h3><p>Computational tours that put the book’s mathematical ideas into practice through detailed implementations.</p><a href="https://www.numerical-tours.com/">Explore the Numerical Tours <span aria-hidden="true">↗</span></a></article>
        <article class="resource-card"><span class="resource-type">Course notes</span><h3>Optimization for ML</h3><p>Focused notes on optimization and its applications in machine learning.</p><ul><li><a href="{{ '/book-sources/aux/optim-ml/OptimML.pdf' | relative_url }}">Optimization for machine learning</a><span class="file-type">PDF</span></li><li><a href="{{ '/book-sources/aux/optim-ml/OptimML-PAC.pdf' | relative_url }}">PAC learning</a><span class="file-type">PDF</span></li><li><a href="{{ '/book-sources/aux/optim-ml/OptimML-DiffusionModels.pdf' | relative_url }}">Diffusion models</a><span class="file-type">PDF</span></li></ul></article>
        <article class="resource-card"><span class="resource-type">Related book</span><h3>Optimal Transport</h3><p>OT4ML: <em>Optimal Transport for Machine Learners</em>. A dedicated book with an interactive edition and teaching resources.</p><a href="https://www.gpeyre.com/ot4ml/">Explore OT4ML <span aria-hidden="true">↗</span></a></article>
      </div>
    </div>
  </section>
</main>
