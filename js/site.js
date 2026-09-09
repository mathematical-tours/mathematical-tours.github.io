/* Navigation and catalogue filtering progressively enhance the static HTML. */
(() => {
  const nav = document.querySelector('.site-nav');
  const toggle = document.querySelector('.nav-toggle');
  const links = document.getElementById('site-navigation');
  const resources = document.querySelector('.nav-resources');
  if (nav && toggle && links) {
    nav.classList.add('nav-enhanced');
    toggle.hidden = false;
    const close = () => {
      toggle.setAttribute('aria-expanded', 'false');
      links.classList.remove('is-open');
      if (resources) resources.open = false;
    };
    toggle.addEventListener('click', () => {
      const open = toggle.getAttribute('aria-expanded') !== 'true';
      toggle.setAttribute('aria-expanded', String(open));
      links.classList.toggle('is-open', open);
    });
    document.addEventListener('keydown', event => {
      if (event.key !== 'Escape') return;
      if (resources && resources.open) {
        resources.open = false;
        resources.querySelector('summary').focus();
      } else if (toggle.getAttribute('aria-expanded') === 'true') {
        close();
        toggle.focus();
      }
    });
    document.addEventListener('click', event => {
      if (resources && !resources.contains(event.target)) resources.open = false;
      if (!nav.contains(event.target)) close();
    });
    links.querySelectorAll('a').forEach(link => link.addEventListener('click', close));
  }

  const controls = document.querySelector('.catalogue-tools');
  if (!controls) return;
  const search = document.getElementById('algorithm-search');
  const topic = document.getElementById('algorithm-topic');
  const rows = [...document.querySelectorAll('[data-algorithm]')];
  const groups = [...document.querySelectorAll('[data-algorithm-group]')];
  const count = document.getElementById('algorithm-count');
  const empty = document.getElementById('no-algorithms');
  const normalize = value => value.normalize('NFD').replace(/[\u0300-\u036f]/g, '').toLowerCase();
  function filter() {
    const terms = normalize(search.value).trim().split(/\s+/).filter(Boolean);
    let matches = 0;
    rows.forEach(row => {
      row.hidden = !terms.every(term => normalize(row.dataset.search).includes(term)) ||
        (topic.value !== '' && row.dataset.topic !== topic.value);
      if (!row.hidden) matches++;
    });
    groups.forEach(group => { group.hidden = !group.querySelector('[data-algorithm]:not([hidden])'); });
    document.querySelectorAll('[data-topic-link]').forEach(link => {
      link.hidden = document.getElementById(link.dataset.topicLink).hidden;
    });
    count.textContent = `${matches} ${matches === 1 ? 'entry' : 'entries'}`;
    empty.hidden = matches !== 0;
  }
  search.addEventListener('input', filter);
  topic.addEventListener('change', filter);
  document.getElementById('reset-algorithms').addEventListener('click', () => {
    search.value = ''; topic.value = ''; filter(); search.focus();
  });
  controls.hidden = false;
  filter();
})();
