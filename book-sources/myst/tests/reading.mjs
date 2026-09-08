/** Check embedded sizing, keyboard access, and reading on a static HTTP server. */
import {chromium} from 'playwright';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const base=(process.env.BOOK_URL||'http://127.0.0.1:8790').replace(/\/$/,'');
const out=path.resolve(root,'../build/interactive-book/browser');fs.mkdirSync(out,{recursive:true});
const browser=await chromium.launch({headless:true,...(process.env.BROWSER_EXECUTABLE?{executablePath:process.env.BROWSER_EXECUTABLE}:{})});
const context=await browser.newContext({viewport:{width:1440,height:1000},reducedMotion:'reduce'});
const issues=[],external=new Set(),checks=[];
await context.route('**/*',route=>{const u=new URL(route.request().url());if(u.protocol.startsWith('http')&&u.origin!==new URL(base).origin){external.add(u.href);return route.abort();}return route.continue();});
const page=await context.newPage();page.on('pageerror',e=>issues.push(String(e)));page.on('response',r=>{if(r.status()>=400)issues.push(`${r.status()}: ${r.url()}`);});
const info=await (await context.request.get(base+'/build-info.json')).json();
async function frameFits(locator){
  await locator.scrollIntoViewIfNeeded();
  await page.waitForFunction(selector=>{
    const iframe=document.querySelector(selector),body=iframe?.contentDocument?.body;
    return body?.dataset.ready&&Math.abs(iframe.getBoundingClientRect().height-body.getBoundingClientRect().height)<=12;
  },await locator.evaluate(e=>'iframe[src="'+e.getAttribute('src')+'"]'),{timeout:15000});
}
try{
  await page.goto(base+'/',{waitUntil:'networkidle'});
  assert.ok((await page.title()).includes('Mathematical Foundations'));
  await page.screenshot({path:path.join(out,'home.png')});checks.push('Homepage');
  await page.goto(base+'/fourier/',{waitUntil:'networkidle'});
  const iframe=page.locator('iframe[src*="fourier--filtering-1d.html"]');await frameFits(iframe);
  const frame=await (await iframe.elementHandle()).contentFrame();
  await frame.locator('.chart-table summary').first().click();await frameFits(iframe);
  assert.equal(await frame.locator('.chart-table[open] table').count(),1);
  await frame.locator('.chart-table summary').first().click();
  await frame.getByRole('tab',{name:'Book figure',exact:true}).click();await frameFits(iframe);
  await frame.getByRole('tab',{name:'Explore',exact:true}).click();await frameFits(iframe);
  await page.screenshot({path:path.join(out,'embedded-desktop.png')});checks.push('Embedded figure, numerical table, and comparison resizing');
  await page.setViewportSize({width:390,height:844});await frameFits(iframe);
  assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth+2),false);
  await page.screenshot({path:path.join(out,'embedded-mobile.png')});checks.push('Embedded mobile sizing');
  await page.goto(base+'/wavelets/',{waitUntil:'networkidle'});
  const contents=page.locator('.myst-outline-collapsible');assert.equal(await contents.getAttribute('aria-expanded'),'false');
  await contents.click();assert.equal(await contents.getAttribute('aria-expanded'),'true');await contents.click();
  await page.getByRole('button',{name:'Search the book'}).click();await page.locator('#search-input').fill('Huffman');
  await page.locator('#search-list a').first().waitFor();assert.ok(await page.locator('#search-list a[href*="shannon-coding"]').count());
  await page.keyboard.press('Escape');checks.push('Mobile contents and accessible search');
  await page.setViewportSize({width:1440,height:1000});await page.goto(base+'/shannon-coding/',{waitUntil:'networkidle'});
  await page.evaluate(()=>document.documentElement.classList.add('dark'));
  await page.screenshot({path:path.join(out,'reading-dark.png')});checks.push('Dark reading theme');
}catch(error){issues.push(String(error));}
await context.close();await browser.close();
const report={book_url:base,content_sha256:info.content_sha256,checks,issues,external_requests:[...external]};
fs.writeFileSync(path.join(out,'reading-report.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report,null,2));
if(issues.length||external.size)process.exitCode=1;
