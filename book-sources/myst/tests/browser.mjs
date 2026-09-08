/** End-to-end review of the built static book, including offline operation. */
import {chromium} from 'playwright';
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import assert from 'node:assert/strict';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const base=(process.env.BOOK_URL||'http://127.0.0.1:8790').replace(/\/$/,'');
const origin=new URL(base).origin;
const out=path.resolve(root,'../build/interactive-book/browser');fs.mkdirSync(out,{recursive:true});
const figures=JSON.parse(fs.readFileSync(path.join(root,'live/figures.json')));
const chapters=JSON.parse(fs.readFileSync(path.join(root,'conversion-report.json'))).chapters;
const executable=process.env.BROWSER_EXECUTABLE;
const browser=await chromium.launch({headless:true,...(executable?{executablePath:executable}:{})});
const context=await browser.newContext({viewport:{width:1280,height:1000},reducedMotion:'reduce'});
const buildInfo=await (await context.request.get(base+'/build-info.json')).json();
const issues=[],results=[],external=new Set();
function saveReport(chaptersChecked=0){const report={book_url:base,content_sha256:buildInfo.content_sha256,figures_checked:results.length,passed:results.filter(r=>r.passed).length,chapters_checked:chaptersChecked,issues,external_requests:[...external],results};fs.writeFileSync(path.join(out,'report.json'),JSON.stringify(report,null,2)+'\n');return report;}
await context.route('**/*',route=>{const u=new URL(route.request().url());if(u.protocol.startsWith('http')&&u.origin!==origin){external.add(u.href);return route.abort();}return route.continue();});
function watch(page,label){page.on('pageerror',e=>issues.push({label,error:String(e)}));page.on('response',r=>{if(r.status()>=400)issues.push({label,status:r.status(),url:r.url()});});}
async function waitPaint(page){await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));}
async function inspectFigure(f){
  const page=await context.newPage();watch(page,f.id);
  try{
    await page.goto(`${base}/live/panels/${f.id}.html`,{waitUntil:'networkidle'});
    await page.waitForFunction(()=>document.body.dataset.ready||document.body.dataset.error,{timeout:15000});
    assert.equal(await page.locator('body').getAttribute('data-error'),null);
    if(f.interaction){
      const initial=await page.locator('body').getAttribute('data-state');const visual=()=>page.locator('#panels').evaluate(p=>p.innerHTML+[...p.querySelectorAll('canvas')].map(c=>c.toDataURL()).join('|'));const before=await visual();
      await page.locator('#controls input,#controls select').first().evaluate(input=>{input.value=input.value==='0'?'1':'0';input.dispatchEvent(new Event('input',{bubbles:true}));});await waitPaint(page);
      assert.notEqual(await page.locator('body').getAttribute('data-state'),initial);
      assert.notEqual(await visual(),before);
      await page.getByRole('button',{name:'Reset all parameters'}).click();await waitPaint(page);assert.equal(await page.locator('body').getAttribute('data-state'),initial);
      await page.getByRole('tab',{name:'Book figure',exact:true}).click();assert.equal(await page.locator('#static').isVisible(),true);
      await page.getByRole('tab',{name:'Explore',exact:true}).click();assert.equal(await page.locator('#explore').isVisible(),true);
      const failed=await page.locator('img:visible').evaluateAll(images=>images.filter(i=>i.complete&&i.naturalWidth===0).map(i=>i.src));assert.deepEqual(failed,[]);
    }else {assert.equal(await page.locator('#static').isVisible(),true);await page.locator('#static img').evaluate(img=>{if(!img.complete||!img.naturalWidth)throw new Error('Static image missing');});}
    results.push({id:f.id,passed:true});
  }catch(error){issues.push({label:f.id,error:String(error)});results.push({id:f.id,passed:false});}
  finally{await page.close();}
}
let next=0;await Promise.all(Array.from({length:4},async()=>{while(next<figures.length){const f=figures[next++];await inspectFigure(f);if(results.length%25===0){saveReport();console.log(`Checked ${results.length}/${figures.length} figure panels`);}}}));
saveReport();
const reader=await context.newPage();watch(reader,'book');
let chaptersChecked=0;
try{
for(const c of chapters){
  try{
  await reader.goto(`${base}/${c.chapter}/`,{waitUntil:'networkidle'});await reader.evaluate(()=>document.fonts.ready);await waitPaint(reader);
  assert.equal(await reader.locator('.katex-error').count(),0,`Math errors in ${c.chapter}`);
  const overflow=await reader.evaluate(()=>document.documentElement.scrollWidth>window.innerWidth+2);assert.equal(overflow,false,`Horizontal overflow in ${c.chapter}`);
  const expected=figures.filter(f=>f.chapter===c.chapter).length;assert.equal(await reader.locator('iframe').count(),expected,`Dropped figure in ${c.chapter}`);
  chaptersChecked++;
  }catch(error){issues.push({label:c.chapter,error:String(error)});}
}
// Real chapter navigation (the MyST router), a mathematical cross-reference,
// and the book search are checked separately from direct static page requests.
await reader.goto(`${base}/shannon-coding/`,{waitUntil:'networkidle'});
await reader.locator('h1').screenshot({path:path.join(out,'chapter-title.png')});
await reader.screenshot({path:path.join(out,'reading-desktop.png')});
await reader.locator('a[href*="fourier"]').first().click();await reader.waitForURL(/fourier/);assert.ok((await reader.title()).includes('Fourier'));
await reader.goto(`${base}/shannon-coding/#eq-proof-shannon`,{waitUntil:'networkidle'});
await reader.locator('#eq-proof-shannon').scrollIntoViewIfNeeded();await reader.screenshot({path:path.join(out,'shannon-proof.png')});
await reader.getByRole('button',{name:/Search/}).first().click();
await reader.locator('#search-input').fill('Shannon');await reader.locator('#search-list a').first().waitFor();
assert.ok(await reader.locator('#search-list a[href*="shannon-coding"]').count());
await reader.locator('#search-input').fill('wavelet');await reader.locator('#search-list a[href*="wavelets"]').first().waitFor();
await reader.locator('#search-list a[href*="wavelets"]').first().click();await reader.waitForURL(/wavelets/);
const samples=['wavelets--wavcoefs-2d','approximation--triangulation-anisoiso','variational-priors--discrete-operators','compressed-sensing--single-pixel','machine-learning--digits-classes','optim-ml-advanced--compariso-sgd'];
for(const id of samples){
  await reader.setViewportSize({width:1000,height:1100});await reader.goto(`${base}/live/panels/${id}.html`,{waitUntil:'networkidle'});await reader.waitForFunction(()=>document.body.dataset.ready==='interactive');await waitPaint(reader);await reader.screenshot({path:path.join(out,id+'.png'),fullPage:true});
  await reader.setViewportSize({width:390,height:844});await waitPaint(reader);assert.equal(await reader.evaluate(()=>document.documentElement.scrollWidth>innerWidth+2),false,`Mobile overflow in ${id}`);await reader.screenshot({path:path.join(out,id+'-mobile.png'),fullPage:true});
}
await reader.goto(`${base}/wavelets/`,{waitUntil:'networkidle'});assert.equal(await reader.evaluate(()=>document.documentElement.scrollWidth>innerWidth+2),false,'Mobile chapter overflow');await reader.screenshot({path:path.join(out,'reading-mobile.png')});
}catch(error){issues.push({label:'reader navigation and layout',error:String(error)});}
await context.close();await browser.close();
const report=saveReport(chaptersChecked);console.log(JSON.stringify({...report,results:undefined},null,2));
if(issues.length||external.size)process.exitCode=1;
