---
book_page: true
social_image: /img/book-basics-cover.webp
social_image_alt: "Cover of the book by Gabriel Peyré"
artwork: paysages-2
layout: book-basics
permalink: /book-basics/
title: "An Introduction to Data Sciences"
description: "Gabriel Peyré's illustrated introduction to information theory, image processing, compressed sensing, optimal transport, and neural networks. Complete books and individual chapters in English and French."
---
<main id="main">
  <section class="hero wrap" aria-labelledby="book-title">
    <div>
      <p class="eyebrow">A visual introduction</p>
      <h1 id="book-title">An Introduction to <em>Data Sciences</em></h1>
      <p class="byline">Gabriel Peyré<span class="affiliation">CNRS &amp; DMA, École Normale Supérieure, PSL University</span></p>
      <p class="hero-intro">How do we compress a message, restore an image, compare distributions, or learn from examples? A short, illustrated journey through the mathematics behind these questions.</p>
      <div class="hero-actions">
        <a class="button primary" href="{{ '/book-basics-sources/IntroImaging.pdf' | relative_url }}" hreflang="en">Full book · English PDF <span aria-hidden="true">↗</span></a>
        <a class="button" href="{{ '/book-basics-sources/IntroImaging-FR.pdf' | relative_url }}" hreflang="fr">Full book · French PDF <span aria-hidden="true">↗</span></a>
      </div>
      <p class="format-note">English: 39 pages · 13.3 MB <span aria-hidden="true">/</span> French: 39 pages · 13.4 MB<span class="basics-edition">Revised September 9, 2026</span></p>
    </div>
    <figure class="cover">
      <a href="{{ '/book-basics-sources/IntroImaging.pdf' | relative_url }}" aria-label="Open the complete English book PDF"><img src="{{ '/img/book-basics-cover.webp' | relative_url }}" alt="Cover of An Introduction to Data Sciences by Gabriel Peyré, with a landscape painting by Louis Peyré" width="894" height="1263" fetchpriority="high"></a>
      <figcaption>Complete editions in English and French</figcaption>
    </figure>
  </section>
  <section class="wrap basics-chapters" id="chapters" aria-labelledby="contents-title">
    <div class="section-heading"><div><p class="eyebrow">Choose a chapter</p><h2 id="contents-title">From information to learning</h2></div><p>Download a chapter in either language, or select an illustration to browse its figures.</p></div>
    <ol class="chapter-list">
      <li>
        <a class="chapter-visual" href="{{ '/book-basics/figures/' | relative_url }}?chapter=1" aria-label="Browse the figures for chapter 1: Claude Shannon and data compression"><img src="{{ '/book-basics-sources/gallery/shared/trees-thumb.webp' | relative_url }}" alt="Prefix codes and binary trees" width="2600" height="660" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">01</span><div><h3 class="chapter-title">Claude Shannon and data compression</h3><p class="chapter-description">Binary codes, entropy, Huffman coding, and the limits of lossless compression.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/book-basics-sources/chapters-pdf/en/01-information.pdf' | relative_url }}" hreflang="en" aria-label="Chapter 1: Claude Shannon and data compression — English PDF">English PDF<span class="pages">8 p.</span></a><a href="{{ '/book-basics-sources/chapters-pdf/fr/01-information.pdf' | relative_url }}" hreflang="fr" aria-label="Chapter 1: Claude Shannon and data compression — French PDF">French PDF<span class="pages">8 p.</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/book-basics/figures/' | relative_url }}?chapter=2" aria-label="Browse the figures for chapter 2: Image processing"><img src="{{ '/book-basics-sources/gallery/en/canaux-coul-thumb.webp' | relative_url }}" alt="Red, green, and blue channels" width="2600" height="900" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">02</span><div><h3 class="chapter-title">Image processing</h3><p class="chapter-description">Pixels, quantization, filtering, edges, colors, and elementary image transformations.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/book-basics-sources/chapters-pdf/en/02-images.pdf' | relative_url }}" hreflang="en" aria-label="Chapter 2: Image processing — English PDF">English PDF<span class="pages">9 p.</span></a><a href="{{ '/book-basics-sources/chapters-pdf/fr/02-images.pdf' | relative_url }}" hreflang="fr" aria-label="Chapter 2: Image processing — French PDF">French PDF<span class="pages">10 p.</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/book-basics/figures/' | relative_url }}?chapter=3" aria-label="Browse the figures for chapter 3: Sparsity, inverse problems and compressed sensing"><img src="{{ '/book-basics-sources/gallery/shared/inpainting-thumb.webp' | relative_url }}" alt="Reconstructing missing pixels" width="2600" height="900" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">03</span><div><h3 class="chapter-title">Sparsity, inverse problems and compressed sensing</h3><p class="chapter-description">Sparse representations, regularization, and recovery from incomplete measurements.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/book-basics-sources/chapters-pdf/en/03-sampling.pdf' | relative_url }}" hreflang="en" aria-label="Chapter 3: Sparsity, inverse problems and compressed sensing — English PDF">English PDF<span class="pages">7 p.</span></a><a href="{{ '/book-basics-sources/chapters-pdf/fr/03-sampling.pdf' | relative_url }}" hreflang="fr" aria-label="Chapter 3: Sparsity, inverse problems and compressed sensing — French PDF">French PDF<span class="pages">8 p.</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/book-basics/figures/' | relative_url }}?chapter=4" aria-label="Browse the figures for chapter 4: Optimal transport"><img src="{{ '/book-basics-sources/gallery/shared/image-eq-thumb.webp' | relative_url }}" alt="Transferring a color palette" width="2600" height="1287" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">04</span><div><h3 class="chapter-title">Optimal transport</h3><p class="chapter-description">Assignments, transport plans, Wasserstein distances, and applications to distributions.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/book-basics-sources/chapters-pdf/en/04-transport.pdf' | relative_url }}" hreflang="en" aria-label="Chapter 4: Optimal transport — English PDF">English PDF<span class="pages">7 p.</span></a><a href="{{ '/book-basics-sources/chapters-pdf/fr/04-transport.pdf' | relative_url }}" hreflang="fr" aria-label="Chapter 4: Optimal transport — French PDF">French PDF<span class="pages">8 p.</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/book-basics/figures/' | relative_url }}?chapter=5" aria-label="Browse the figures for chapter 5: Neural networks"><img src="{{ '/book-basics-sources/gallery/en/discriminative-thumb.webp' | relative_url }}" alt="An image classifier" width="2600" height="1378" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">05</span><div><h3 class="chapter-title">Neural networks</h3><p class="chapter-description">Classification, backpropagation, generative models, and connections with transport.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/book-basics-sources/chapters-pdf/en/05-neural-networks.pdf' | relative_url }}" hreflang="en" aria-label="Chapter 5: Neural networks — English PDF">English PDF<span class="pages">7 p.</span></a><a href="{{ '/book-basics-sources/chapters-pdf/fr/05-neural-networks.pdf' | relative_url }}" hreflang="fr" aria-label="Chapter 5: Neural networks — French PDF">French PDF<span class="pages">7 p.</span></a></div>
      </li>
    </ol>
  </section>
  <section class="about wrap" id="about" aria-labelledby="about-title">
    <h2 id="about-title">A mathematical starting point</h2>
    <div>
      <p>The first two chapters start from elementary algebra, logarithms, and basic probability. The later chapters introduce more advanced ideas using vectors, matrices, derivatives, and probability distributions. Concrete examples motivate the definitions and theorems.</p>
      <p>Both editions include the five chapters, including discriminative and generative neural networks. This revision aligns the two texts, refines the mathematical statements, and updates the typography.</p>
    </div>
  </section>
  <section class="resources" id="resources" aria-labelledby="resources-title"><div class="wrap">
    <h2 id="resources-title">Figures and further reading</h2>
    <div class="resource-grid resource-grid-two">
      <article class="resource-card"><span class="resource-type">Visual companion</span><h3>Explore all 42 figures</h3><p>Browse by chapter or search by keyword, in English or French. Enlarge each illustration and download its PDF or JPG, with the complete panels and labels.</p><a href="{{ '/book-basics/figures/' | relative_url }}">Browse the figure gallery <span aria-hidden="true">↗</span></a></article>
      <article class="resource-card"><span class="resource-type">Continue reading</span><h3>Mathematical Foundations of Data Sciences</h3><p>A more extensive treatment of Fourier analysis, wavelets, approximation, inverse problems, and optimization, developing the ideas introduced here.</p><a href="{{ '/book/' | relative_url }}">Explore the advanced book <span aria-hidden="true">↗</span></a></article>
    </div>
  </div></section>
</main>
