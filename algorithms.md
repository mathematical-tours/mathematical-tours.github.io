---
layout: resource
permalink: /algorithms/
title: "Algorithms"
description: "A searchable directory of algorithms and mathematical models, grouped by subject: coding, graphs, matrices, geometry, signals, learning, probability, and optimization."
intro: "A directory of methods that shape computational mathematics. Browse by subject or search for an algorithm, then follow its reference to learn more."
resource_page: true
artwork: paysages-paris-1
---
{% assign total = 0 %}{% for group in site.data.algorithms %}{% assign total = total | plus: group.entries.size %}{% endfor %}
<div class="page-meta"><span><strong>{{ total }} entries</strong> to explore</span><span><strong>{{ site.data.algorithms.size }} subjects</strong></span><span>Algorithms, models, and computational ideas</span></div>
<div class="catalogue-tools" hidden>
  <div class="field"><label for="algorithm-search">Find a method</label><input type="search" id="algorithm-search" placeholder="Try Fourier, shortest path, or sampling" autocomplete="off" aria-controls="algorithm-catalogue"></div>
  <div class="field"><label for="algorithm-topic">Subject</label><select id="algorithm-topic" aria-controls="algorithm-catalogue"><option value="">All subjects</option>{% for group in site.data.algorithms %}<option value="{{ group.id }}">{{ group.title }}</option>{% endfor %}</select></div>
  <p class="result-count" id="algorithm-count" role="status" aria-live="polite">{{ total }} entries</p>
</div>
<div class="directory-layout">
  <nav class="section-nav" aria-label="Algorithm subjects"><p>Browse by subject</p>{% for group in site.data.algorithms %}<a href="#{{ group.id }}" data-topic-link="{{ group.id }}">{{ group.title }}</a>{% endfor %}</nav>
  <div id="algorithm-catalogue">
    {% for group in site.data.algorithms %}
    <section class="directory-section" id="{{ group.id }}" data-algorithm-group aria-labelledby="title-{{ group.id }}">
      <h2 id="title-{{ group.id }}">{{ group.title }}</h2><p>{{ group.description }}</p>
      <ul class="algorithm-list">{% for item in group.entries %}<li data-algorithm data-topic="{{ group.id }}" data-search="{{ item.name | escape }} {{ group.title | escape }} {{ group.description | escape }}"><a href="{{ item.url }}"><span>{{ item.name | escape }}</span><span aria-hidden="true">↗</span></a></li>{% endfor %}</ul>
    </section>
    {% endfor %}
    <p class="empty-results" id="no-algorithms" hidden>No methods match your search. <button type="button" id="reset-algorithms">Show all entries</button></p>
  </div>
</div>
<p class="editorial-note">This collection includes both algorithms and mathematical models that can be solved in several ways. The subjects overlap: a method is listed in one group to make browsing easier. References provide definitions and historical context; the directory is not a ranking or an exhaustive history. Suggestions and corrections are welcome through <a href="https://github.com/mathematical-tours/mathematical-tours.github.io/issues">the website repository</a>.</p>
