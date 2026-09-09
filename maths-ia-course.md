---
layout: resource
permalink: /maths-ia-course/
title: "Mathematics of Machine Learning"
description: "Ten two-hour lectures on optimization, neural networks, generative models, and optimal transport, with transcripts, notebooks, slides, and reading lists."
intro: "A course in 10 lectures, from optimization to generative models and optimal transport. Follow the mathematical ideas, then explore them through notebooks and worked examples."
resource_page: true
banner: /img/boat-waves.jpg
---
<div class="page-meta"><span><strong>10 lectures</strong> · 2 hours each</span><span><strong>9 PDF transcripts</strong></span><span>Notebooks · Slides · Further reading</span></div>
<div class="course-introduction"><p><strong>Before you begin.</strong> Familiarity with linear algebra, multivariable calculus, and basic probability will help. The lectures emphasize the main concepts and methods; the linked books provide fuller proofs and background.</p><p>Each lecture gathers its topics, practical materials, and suggested reading. Transcripts are available for lectures 1–9. Lecture 10 is supported by the <a href="https://www.gpeyre.com/ot4ml/">OT4ML book</a>, slides, and notebooks.</p></div>
<div class="directory-layout">
  <nav class="section-nav" aria-label="Lecture contents"><p>The ten lectures</p>{% for lecture in site.data.lectures %}<a href="#lecture-{{ lecture.number }}">{% if lecture.number < 10 %}0{% endif %}{{ lecture.number }} · {{ lecture.title }}</a>{% endfor %}</nav>
  <div>{% for lecture in site.data.lectures %}
    <section class="lecture" id="lecture-{{ lecture.number }}" aria-labelledby="lecture-title-{{ lecture.number }}">
      <p class="eyebrow">Lecture {% if lecture.number < 10 %}0{% endif %}{{ lecture.number }}</p><h2 id="lecture-title-{{ lecture.number }}">{{ lecture.title }}</h2>
      {% if lecture.transcript %}<a class="button" href="{{ lecture.transcript | relative_url }}" aria-label="Lecture {{ lecture.number }}: {{ lecture.title }} — PDF transcript">Read the transcript <span class="file-type">PDF</span><span aria-hidden="true">↗</span></a>{% else %}<a class="button" href="https://www.gpeyre.com/ot4ml/">Read OT4ML <span aria-hidden="true">↗</span></a>{% endif %}
      <div class="lecture-columns"><div><h3>Topics</h3>{{ lecture.content | markdownify }}</div><div><h3>Notebooks &amp; materials</h3>{{ lecture.materials | markdownify }}</div></div>
      {% if lecture.bibliography %}<div class="lecture-bibliography"><h3>Further reading</h3>{{ lecture.bibliography | markdownify }}</div>{% endif %}
    </section>
  {% endfor %}</div>
</div>
