---
book_page: true
social_image: /img/daft-cover.webp
social_image_alt: "Cover of the book by Gabriel Peyré"
artwork: natures-mortes-2
layout: daft
permalink: /daft/
title: "The Discrete Algebra of the Fourier Transform"
description: "Gabriel Peyré's DAFT book: an adapted English translation of L'algèbre discrète de la transformée de Fourier (Ellipses, 2004), with Python programs, notebooks, and a searchable figure gallery."
---
<main id="main">
  <section class="hero wrap" aria-labelledby="book-title">
    <div>
      <p class="eyebrow">The DAFT book · Revised English edition</p>
      <h1 id="book-title">The Discrete Algebra of the <em>Fourier Transform</em></h1>
      <p class="byline">Gabriel Peyré<span class="affiliation">CNRS &amp; DMA, École Normale Supérieure, PSL University</span></p>
      <p class="hero-intro">The algebra behind Fourier analysis, from finite groups and fast algorithms to error-correcting codes and representations. Mathematical ideas developed through proofs, exercises, and working Python examples.</p>
      <div class="hero-actions">
        <a class="button primary" href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}">Read the book <span aria-hidden="true">↗</span></a>
        <a class="button" href="{{ '/daft/figures/' | relative_url }}">Explore the figures <span aria-hidden="true">↗</span></a>
      </div>
      <p class="format-note">PDF · 206 pages · 4.2 MB · September 2026</p>
    </div>
    <figure class="cover">
      <a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}" aria-label="Open the DAFT book PDF"><img src="{{ '/img/daft-cover.webp' | relative_url }}" alt="The Discrete Algebra of the Fourier Transform — Gabriel Peyré" width="893" height="1263" fetchpriority="high"></a>
      <figcaption>Adapted and translated from the author's French book<span><em>L'algèbre discrète de la transformée de Fourier</em> · Ellipses, 2004</span></figcaption>
    </figure>
  </section>
  <section class="wrap daft-chapters" id="chapters" aria-labelledby="contents-title">
    <div class="section-heading"><div><p class="eyebrow">Explore the book</p><h2 id="contents-title">From characters to computation</h2></div><p>Read a chapter in the PDF, or select an illustration to explore its mathematics. Each chapter ends with its exercise solutions; two Python appendices complete the book.</p></div>
    <ol class="chapter-list">
      <li>
        <a class="chapter-visual" href="{{ '/daft/figures/' | relative_url }}?chapter=1#fig-cyclic-characters" aria-label="Browse the figures for chapter 1: Fourier analysis on finite groups"><img src="{{ '/img/daft-figures/chapter-finite-groups/cyclic-characters-thumb.webp' | relative_url }}" alt="Characters of a cyclic group" width="2001" height="855" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">01</span><div><h3 class="chapter-title">Fourier analysis on finite groups</h3><p class="chapter-description">Characters, duality, orthogonality, convolution, and random walks.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page=5" aria-label="Read chapter 1: Fourier analysis on finite groups">Read chapter <span aria-hidden="true">↗</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/daft/figures/' | relative_url }}?chapter=2#fig-functions-walsh-2d" aria-label="Browse the figures for chapter 2: Applications of Fourier duality"><img src="{{ '/img/daft-figures/chapter-duality/walsh-basis-two-dimensions-thumb.webp' | relative_url }}" alt="Two-dimensional Walsh functions" width="1600" height="1600" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">02</span><div><h3 class="chapter-title">Applications of Fourier duality</h3><p class="chapter-description">Gauss sums, quadratic reciprocity, Walsh transforms, and Poisson summation.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page=25" aria-label="Read chapter 2: Applications of Fourier duality">Read chapter <span aria-hidden="true">↗</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/daft/figures/' | relative_url }}?chapter=3#fig-fft-butterfly-8points" aria-label="Browse the figures for chapter 3: The discrete Fourier transform"><img src="{{ '/img/daft-figures/chapter-discrete-fourier-transform/eight-point-fft-thumb.webp' | relative_url }}" alt="An eight-point FFT" width="2000" height="976" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">03</span><div><h3 class="chapter-title">The discrete Fourier transform</h3><p class="chapter-description">FFT factorizations, multidimensional transforms, interpolation, and spectral structure.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page=54" aria-label="Read chapter 3: The discrete Fourier transform">Read chapter <span aria-hidden="true">↗</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/daft/figures/' | relative_url }}?chapter=4#fig-fourier-descriptors" aria-label="Browse the figures for chapter 4: Applications of the discrete Fourier transform"><img src="{{ '/img/daft-figures/chapter-fourier-applications/fourier-descriptors-thumb.webp' | relative_url }}" alt="Fourier shape descriptors" width="2000" height="898" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">04</span><div><h3 class="chapter-title">Applications of the discrete Fourier transform</h3><p class="chapter-description">Filtering, geometric descriptors, differential equations, and fast products.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page=82" aria-label="Read chapter 4: Applications of the discrete Fourier transform">Read chapter <span aria-hidden="true">↗</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/daft/figures/' | relative_url }}?chapter=5#fig-location-poles-notch-filter" aria-label="Browse the figures for chapter 5: Transforms related to the DFT"><img src="{{ '/img/daft-figures/chapter-related-transforms/notch-poles-and-zeros-thumb.webp' | relative_url }}" alt="Poles and zeros of a notch filter" width="1645" height="1600" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">05</span><div><h3 class="chapter-title">Transforms related to the DFT</h3><p class="chapter-description">Hartley transforms, the z-transform, chirp algorithms, and scaled frequencies.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page=113" aria-label="Read chapter 5: Transforms related to the DFT">Read chapter <span aria-hidden="true">↗</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/daft/figures/' | relative_url }}?chapter=6#fig-visual-dist-hamming" aria-label="Browse the figures for chapter 6: Fourier transforms over finite fields and rings"><img src="{{ '/img/daft-figures/chapter-finite-fields-and-codes/hamming-cube-thumb.webp' | relative_url }}" alt="Hamming distance on the cube" width="1913" height="1600" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">06</span><div><h3 class="chapter-title">Fourier transforms over finite fields and rings</h3><p class="chapter-description">Exact transforms, Hamming and BCH codes, decoding, and Boolean functions.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page=135" aria-label="Read chapter 6: Fourier transforms over finite fields and rings">Read chapter <span aria-hidden="true">↗</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/daft/figures/' | relative_url }}?chapter=7#wavelets-finite-body" aria-label="Browse the figures for chapter 7: Linear representations of finite groups"><img src="{{ '/img/daft-figures/chapter-linear-representations/finite-field-wavelets-thumb.webp' | relative_url }}" alt="Wavelets over a finite field" width="1890" height="1600" loading="lazy" decoding="async"><span>View figures <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">07</span><div><h3 class="chapter-title">Linear representations of finite groups</h3><p class="chapter-description">Irreducibility, characters, matrix-valued Fourier analysis, and finite-field wavelets.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page=160" aria-label="Read chapter 7: Linear representations of finite groups">Read chapter <span aria-hidden="true">↗</span></a></div>
      </li>
      <li>
        <a class="chapter-visual" href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page=181" aria-label="Read about cube symmetries in chapter 8"><img src="{{ '/img/previews/daft-cube-symmetries.svg' | relative_url }}" alt="Four body diagonals permuted by the 24 rotations of a cube" width="440" height="240" loading="lazy" decoding="async"><span>Cube symmetries <span aria-hidden="true">↗</span></span></a>
        <div class="chapter-copy"><span class="chapter-number">08</span><div><h3 class="chapter-title">Applications of linear representations</h3><p class="chapter-description">Character tables, cube symmetries, ranking data, quaternions, and invariant polynomials.</p></div></div>
        <div class="chapter-actions"><a href="{{ '/daft-sources/DAFT-EN.pdf' | relative_url }}#page=179" aria-label="Read chapter 8: Applications of linear representations">Read chapter <span aria-hidden="true">↗</span></a></div>
      </li>
    </ol>
  </section>
  <section class="about wrap" id="edition" aria-labelledby="edition-title">
    <h2 id="edition-title">An English adaptation of the French original</h2>
    <div><p>This book is an adapted English translation of Gabriel Peyré's <a href="https://www.editions-ellipses.fr/fr/accueil/12783-23434-l-algebre-discrete-de-la-transformee-de-fourier-niveau-m1-9782729818678.html" lang="fr"><em>L'algèbre discrète de la transformée de Fourier</em></a>, published by <strong>Ellipses in 2004</strong>, ISBN 9782729818678.</p><p>The September 2026 edition develops the proofs and exposition with 68 figures and Python examples throughout. It assumes undergraduate linear algebra and basic familiarity with groups and complex numbers.</p><p>The previous English edition and its original programs remain available in the <a href="https://github.com/mathematical-tours/mathematical-tours.github.io/blob/master/daft-sources/legacy/pre-python-2026-09-08/DAFT-EN.pdf">legacy PDF</a> and <a href="https://github.com/mathematical-tours/mathematical-tours.github.io/tree/master/daft-sources/legacy">archived sources</a>.</p></div>
  </section>
  <section class="resources" id="resources" aria-labelledby="resources-title"><div class="wrap">
    <h2 id="resources-title">Figures and Python</h2>
    <div class="resource-grid">
      <article class="resource-card"><span class="resource-type">Visual companion</span><h3>Explore all 68 figures</h3><p>Search by topic, chapter, or figure number. Open a full-size preview, download the vector PDF, and follow the Python notebook or TikZ source.</p><a href="{{ '/daft/figures/' | relative_url }}">Browse the figure gallery <span aria-hidden="true">↗</span></a></article>
      <article class="resource-card"><span class="resource-type">Companion code</span><h3>Python examples</h3><p>The book's numerical and symbolic programs, with shared helpers and installation instructions. Use them to explore transforms, finite groups, interpolation, and error-correcting codes.</p><a href="{{ '/daft-sources/downloads/daft-python.zip' | relative_url }}" download>Download the Python examples <span aria-hidden="true">↓</span></a></article>
      <article class="resource-card"><span class="resource-type">Reproduce the figures</span><h3>Code and data</h3><p>All figure generators in one download, including the plotting helpers, TikZ diagrams, and flower photograph used in the image examples.</p><a href="{{ site.data.daft_figures.archive | relative_url }}" download>Download figure code and data <span aria-hidden="true">↓</span></a></article>
    </div>
  </div></section>
</main>
