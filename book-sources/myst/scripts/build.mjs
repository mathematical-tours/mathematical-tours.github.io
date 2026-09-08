import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import path from 'node:path';
import fs from 'node:fs';
import {createHash} from 'node:crypto';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const base=(process.env.BASE_URL||'').replace(/\/$/,'');
const site=(process.env.SITE_URL||'https://mathematical-tours.github.io').replace(/\/$/,'');
const publicURL=site+base;
// MyST retains hashed downloads between builds. Recreate only the generated
// site/export trees so a publication never contains obsolete chapter sources.
for(const name of ['html','site'])fs.rmSync(path.join(root,'_build',name),{recursive:true,force:true});
const args=['build','--html'];if(process.argv.includes('--strict'))args.push('--strict');
const result=spawnSync(path.join(root,'node_modules/.bin/myst'),args,{cwd:root,stdio:'inherit',env:{...process.env,BASE_URL:base,HOST:'127.0.0.1'}});
if(result.status!==0)process.exit(result.status||1);
function walk(p){return fs.readdirSync(p,{withFileTypes:true}).flatMap(e=>e.isDirectory()?walk(path.join(p,e.name)):[path.join(p,e.name)]);}
// The upstream theme references three CDN stylesheets even without notebooks.
// Host its matching math CSS/fonts locally; the two widget-only stylesheets are
// unused by our native panels. Patch both HTML and the hydrated client bundle.
const vendor=path.join(root,'_build/html/vendor/katex');fs.mkdirSync(vendor,{recursive:true});
const katexPackage=path.join(root,'node_modules/katex');
const katexCSS=fs.readFileSync(path.join(katexPackage,'dist/katex.min.css'));
fs.writeFileSync(path.join(vendor,'katex.min.css'),katexCSS);
fs.cpSync(path.join(katexPackage,'dist/fonts'),path.join(vendor,'fonts'),{recursive:true});
fs.copyFileSync(path.join(katexPackage,'LICENSE'),path.join(vendor,'LICENSE'));
const integrity='sha384-'+createHash('sha384').update(katexCSS).digest('base64');
const previousIntegrity='sha384-MlJdn/WNKDGXveldHDdyRP1R4CTHr3FeuDNfhsLPYrq2t0UBkUdK2jyTnXPEK1NQ';
const externalMath='https://cdn.jsdelivr.net/npm/katex@0.15.2/dist/katex.min.css';
const unusedStyles=['https://cdn.jsdelivr.net/npm/jupyter-matplotlib@0.11.3/css/mpl_widget.css','https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.css'];
// MyST caches hashed copies of every download and also copies static_files into
// build/. Share one authoritative PDF tree and retain chapter-to-book links.
const output=path.join(root,'_build/html'),pdfRoot=path.join(output,'assets/pdfs');
const pdfByHash=new Map(walk(pdfRoot).filter(p=>p.endsWith('.pdf')).map(p=>[createHash('sha256').update(fs.readFileSync(p)).digest('hex'),path.relative(pdfRoot,p).split(path.sep).join('/')]));
const duplicatePDFs=new Map();
for(const name of fs.readdirSync(path.join(output,'build'))){
  if(!name.endsWith('.pdf'))continue;
  const hash=createHash('sha256').update(fs.readFileSync(path.join(output,'build',name))).digest('hex');
  if(pdfByHash.has(hash))duplicatePDFs.set(name,pdfByHash.get(hash));
}
const referencedCopies=new Set();
const markdownDownloads=new Set(fs.readdirSync(path.join(output,'build')).filter(name=>name.endsWith('.md')));
let count=0;for(const p of walk(path.join(root,'_build/html'))){
  if(p.includes('/live/')||!(/\.(html|json|js)$/.test(p)))continue;
  let text=fs.readFileSync(p,'utf8');
  text=text.replaceAll('"domain":"http://localhost:3000"',`"domain":"${publicURL}"`);
  // Offer the publication PDFs in the theme's download menu. Generated MyST
  // sources have YAML front matter, which the parent Jekyll site would render
  // as unrelated HTML pages instead of preserving them as source downloads.
  text=text.replace(/"exports":\[\{"format":"md","filename":"([^"]+)\.md","url":"[^"]+\.md"\}\]/g,(_,slug)=>{
    const file=['index','references'].includes(slug)?'FundationsDataScience.pdf':`chapters-pdf/${slug}.pdf`;
    if(!fs.existsSync(path.join(pdfRoot,file)))throw new Error(`Missing publication download: ${file}`);
    return '"exports":'+JSON.stringify([{format:'pdf',filename:path.basename(file),url:`${base}/assets/pdfs/${file}`}]);
  });
  text=text.replace(/build\/([\w.-]+\.pdf)/g,(match,name)=>duplicatePDFs.has(name)?'assets/pdfs/'+duplicatePDFs.get(name):match);
  text=text.replaceAll(externalMath,`${base}/vendor/katex/katex.min.css`).replaceAll(previousIntegrity,integrity);
  for(const url of unusedStyles){
    text=text.replaceAll(`<link rel="stylesheet" href="${url}"/>`,'');
    text=text.replaceAll(`{rel:"stylesheet",href:"${url}"},`,'').replaceAll(`{rel:"stylesheet",href:"${url}"}`,'');
  }
  // Rooted iframe paths also work after client-side navigation to a URL without
  // a trailing slash, including deployment below an arbitrary BASE_URL.
  text=text.replace(/\.\.\/live\/panels\//g,`${base}/live/panels/`);
  if(p.endsWith('.html')){
    const page=path.relative(output,path.dirname(p)).split(path.sep).join('/');
    text=text.replace('</head>',`<link rel="canonical" href="${publicURL}/${page?page+'/':''}"/></head>`);
    text=text.replace(/(<meta property="og:image" content=")\//g,`$1${site}/`);
    text=text.replace(/<iframe([^>]*src="[^"]*live\/panels\/[^>]+)>/g,'<iframe loading="lazy"$1>');
    const tag=`<script src="${base}/book-frame.js" defer data-foundations-frames></script>`;
    if(!text.includes('data-foundations-frames'))text=text.replace('</body>',tag+'</body>');
    count++;
  }
  for(const name of ['live','assets'])if(text.includes(`build/${name}/`))referencedCopies.add(name);
  for(const name of markdownDownloads)if(text.includes(name))throw new Error(`Unconverted source download: ${name}`);
  fs.writeFileSync(p,text);
}
for(const name of markdownDownloads)fs.unlinkSync(path.join(output,'build',name));
for(const name of duplicatePDFs.keys())fs.unlinkSync(path.join(output,'build',name));
for(const name of ['live','assets'])if(!referencedCopies.has(name))fs.rmSync(path.join(output,'build',name),{recursive:true,force:true});
for(const name of ['sitemap.xml','robots.txt']){
  const p=path.join(output,name);
  if(fs.existsSync(p))fs.writeFileSync(p,fs.readFileSync(p,'utf8').replaceAll('http://localhost:3000',publicURL));
}
const fingerprint=createHash('sha256');
fingerprint.update(JSON.stringify({base,site}));
const inputs=[...walk(path.join(root,'chapters')),...walk(path.join(root,'live')).filter(p=>!p.includes('/data/images/')), ...['index.md','references.md','myst.yml','book-frame.js','book-theme.css','assets/figure-hashes.json','package-lock.json','scripts/build.mjs'].map(p=>path.join(root,p))].sort();
for(const p of inputs){fingerprint.update(path.relative(root,p));fingerprint.update(fs.readFileSync(p));}
fs.writeFileSync(path.join(root,'_build/html/build-info.json'),JSON.stringify({base_url:base,public_url:publicURL+'/',book:'Mathematical Foundations of Data Sciences',content_sha256:fingerprint.digest('hex')})+'\n');
console.log(`Installed responsive figures on ${count} book pages.`);
