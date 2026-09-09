/* The complete gallery remains available when JavaScript is disabled. */
(() => {
  const cards = [...document.querySelectorAll('[data-figure]')];
  if (!cards.length) return;
  const search = document.getElementById('figure-search');
  const chapter = document.getElementById('figure-chapter');
  const type = document.getElementById('figure-type');
  const count = document.getElementById('figure-count');
  const reset = document.getElementById('reset-figures');
  const normalize = value => value.normalize('NFD').replace(/[\u0300-\u036f]/g, '').toLowerCase();
  const index = new Map(cards.map(card => [card, normalize(card.dataset.search)]));
  const numberQuery = /^(?:fig(?:ure)?\.?\s*)?(\d+\.\d+)\.?$/i;

  function filter(updateURL = true) {
    const query = normalize(search.value).trim();
    const number = query.match(numberQuery);
    const terms = query.split(/\s+/).filter(Boolean);
    let matches = 0;
    for (const card of cards) {
      const textMatch = number ? card.dataset.number === number[1] : terms.every(term => index.get(card).includes(term));
      card.hidden = !textMatch || (chapter.value !== '' && card.dataset.chapter !== chapter.value) ||
        (type.value !== '' && card.dataset.type !== type.value);
      if (!card.hidden) matches++;
    }
    count.textContent = `${matches} of ${cards.length} figures`;
    reset.hidden = !query && !chapter.value && !type.value;
    document.getElementById('no-figures').hidden = matches !== 0;
    if (updateURL) {
      const url = new URL(location.href);
      for (const [key, value] of [['q', search.value.trim()], ['chapter', chapter.value], ['code', type.value]]) {
        if (value) url.searchParams.set(key, value); else url.searchParams.delete(key);
      }
      history.replaceState(null, '', url);
    }
  }

  function restore() {
    const params = new URLSearchParams(location.search);
    search.value = params.get('q') || '';
    chapter.value = params.get('chapter') || '';
    type.value = params.get('code') || '';
    if (chapter.selectedIndex < 0) chapter.value = '';
    if (type.selectedIndex < 0) type.value = '';
    filter(false);
  }
  function clear() {
    search.value = ''; chapter.value = ''; type.value = '';
    filter(); search.focus();
  }
  search.addEventListener('input', () => filter());
  chapter.addEventListener('change', () => filter());
  type.addEventListener('change', () => filter());
  reset.addEventListener('click', clear);
  document.getElementById('reset-empty-figures').addEventListener('click', clear);
  window.addEventListener('popstate', restore);
  document.querySelector('.figure-tools').hidden = false;
  restore();

  const dialog = document.getElementById('figure-viewer');
  if (typeof dialog.showModal !== 'function') return;
  const image = document.getElementById('viewer-image');
  const previous = document.getElementById('previous-figure');
  const next = document.getElementById('next-figure');
  let selection = [], position = 0, opener;
  function show(index) {
    position = index;
    const card = selection[position];
    document.getElementById('viewer-title').textContent = `Figure ${card.dataset.number} · ${card.querySelector('h2').textContent}`;
    document.getElementById('viewer-image-error').hidden = true;
    image.hidden = false;
    image.alt = card.querySelector('img').alt;
    image.src = card.querySelector('[data-enlarge]').dataset.preview;
    document.getElementById('viewer-links').replaceChildren(...[...card.querySelector('.figure-links').children].map(link => link.cloneNode(true)));
    document.getElementById('viewer-position').textContent = `${position + 1} / ${selection.length}`;
    previous.disabled = position === 0;
    next.disabled = position === selection.length - 1;
  }
  image.addEventListener('error', () => {
    image.hidden = true;
    document.getElementById('viewer-image-error').hidden = false;
  });
  document.getElementById('figure-grid').addEventListener('click', event => {
    const link = event.target.closest('[data-enlarge]');
    if (!link || event.ctrlKey || event.metaKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    opener = link;
    selection = cards.filter(card => !card.hidden);
    show(selection.indexOf(link.closest('[data-figure]')));
    dialog.showModal();
    document.body.classList.add('figure-viewer-open');
    dialog.querySelector('.viewer-close').focus();
  });
  previous.addEventListener('click', () => { if (position > 0) show(position - 1); });
  next.addEventListener('click', () => { if (position + 1 < selection.length) show(position + 1); });
  dialog.addEventListener('keydown', event => {
    if (event.key === 'ArrowLeft' && position > 0) { event.preventDefault(); show(position - 1); }
    if (event.key === 'ArrowRight' && position + 1 < selection.length) { event.preventDefault(); show(position + 1); }
  });
  dialog.querySelector('.viewer-close').addEventListener('click', () => dialog.close());
  dialog.addEventListener('click', event => {
    const rect = dialog.getBoundingClientRect();
    if (event.target === dialog && (event.clientX < rect.left || event.clientX > rect.right || event.clientY < rect.top || event.clientY > rect.bottom)) dialog.close();
  });
  dialog.addEventListener('close', () => {
    document.body.classList.remove('figure-viewer-open');
    opener?.focus();
  });
})();
