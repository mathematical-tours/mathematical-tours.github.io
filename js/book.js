/* The chapter catalogue remains fully readable when JavaScript is unavailable. */
const tools=document.querySelector('.chapter-tools');
if(tools){
  const query=document.getElementById('chapter-search'),topic=document.getElementById('chapter-topic');
  const rows=[...document.querySelectorAll('.chapter-list li')],groups=[...document.querySelectorAll('.chapter-group')];
  const clean=text=>text.normalize('NFD').replace(/[\u0300-\u036f]/g,'').toLowerCase();
  function filter(){
    const terms=clean(query.value).trim().split(/\s+/).filter(Boolean);let count=0;
    for(const row of rows){
      const text=clean(row.dataset.search);
      row.hidden=!(terms.every(term=>text.includes(term))&&(!topic.value||row.dataset.topic===topic.value));
      if(!row.hidden)count++;
    }
    for(const group of groups)group.hidden=![...group.querySelectorAll('.chapter-list li')].some(row=>!row.hidden);
    document.getElementById('chapter-count').textContent=`${count} ${count===1?'chapter':'chapters'}`;
    document.getElementById('no-chapters').hidden=count!==0;
  }
  query.addEventListener('input',filter);topic.addEventListener('change',filter);
  document.getElementById('reset-chapters').addEventListener('click',()=>{query.value='';topic.value='';filter();query.focus();});
  tools.hidden=false;filter();
}
