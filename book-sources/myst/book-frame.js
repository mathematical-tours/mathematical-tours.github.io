// MyST wraps embedded iframes in an aspect-ratio container. Give the real
// experiment its measured height, including after client-side chapter changes.
const prepared = new WeakSet();
function prepareFrames(){
  for(const button of document.querySelectorAll('button.myst-search-bar'))button.setAttribute('aria-label','Search the book');
  for(const frame of document.querySelectorAll('.book-live-frame iframe')){
    if(prepared.has(frame))continue;prepared.add(frame);frame.loading='lazy';
    const container=frame.parentElement;container.style.paddingBottom='0';container.style.height='560px';container.style.width='100%';
  }
}
window.addEventListener('message',event=>{
  if(event.origin!==location.origin||event.data?.type!=='foundations-figure-height')return;
  const height=Number(event.data.height);if(!Number.isFinite(height)||height<50||height>12000)return;
  for(const frame of document.querySelectorAll('.book-live-frame iframe')){
    if(frame.contentWindow===event.source){frame.parentElement.style.height=`${height+4}px`;frame.parentElement.style.paddingBottom='0';break;}
  }
});
new MutationObserver(prepareFrames).observe(document.documentElement,{childList:true,subtree:true});prepareFrames();
// Start the native section outline collapsed on phones. Readers can open it
// normally; a remembered button is never automatically closed a second time.
const compactOutlines=new WeakSet();
function prepareMobileOutline(){
  if(!matchMedia('(max-width: 1023px)').matches)return;
  for(const button of document.querySelectorAll('.myst-outline-collapsible')){
    if(compactOutlines.has(button))continue;compactOutlines.add(button);
    if(button.getAttribute('aria-expanded')==='true')requestAnimationFrame(()=>button.click());
  }
}
new MutationObserver(prepareMobileOutline).observe(document.documentElement,{childList:true,subtree:true});
window.addEventListener('resize',prepareMobileOutline);prepareMobileOutline();
