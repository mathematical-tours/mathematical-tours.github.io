---
layout: resource
permalink: /slides/
title: "Course slides"
description: "Seventeen slide decks for teaching signal and image processing, sparse recovery, geometry, and mesh processing."
intro: "Two collections for teaching and review. Follow a course in order, or open the deck for the topic you need."
resource_page: true
banner: /img/hiroshige-whirpool.jpg
---
<div class="page-meta"><span><strong>17 slide decks</strong></span><span><strong>2 course collections</strong></span><span>View the slides on Speaker Deck</span></div>
<nav class="jump-links" aria-label="Slide collections"><a class="button" href="#signals">Signal &amp; image processing ↓</a><a class="button" href="#geometry">Geometry &amp; meshes ↓</a></nav>
{% for group in site.data.slides %}
<section class="slide-section" id="{{ group.id }}" aria-labelledby="title-{{ group.id }}">
  <div class="section-heading"><h2 id="title-{{ group.id }}">{{ group.title }}</h2><p>{{ group.intro }}</p></div>
  <ol class="slide-grid">{% for deck in group.decks %}<li class="slide-card"><span class="slide-number" aria-hidden="true">{% if forloop.index < 10 %}0{% endif %}{{ forloop.index }}</span><div><h3><a href="{{ deck.url }}">{{ deck.title }}</a></h3><p>{{ deck.description }}</p><a class="text-link" href="{{ deck.url }}" aria-label="Open slides: {{ deck.title | escape }}">Open the slides <span aria-hidden="true">↗</span></a></div></li>{% endfor %}</ol>
</section>
{% endfor %}
<div class="editorial-note">For the mathematical development behind the signal-processing course, read <a href="{{ '/book/' | relative_url }}">Mathematical Foundations of Data Sciences</a>. For neural networks, generative models, and optimal transport, follow the <a href="{{ '/maths-ia-course/' | relative_url }}">ten-lecture machine-learning course</a>.</div>
