import {drawPlot,colors,format,plotColor} from './plot.mjs';
import {getGuide,controlHelp,displayValue} from './guidance.mjs';
const $=id=>document.getElementById(id);
const config=JSON.parse($('figure-config').textContent);
const el=(tag,cls,text)=>{const x=document.createElement(tag);if(cls)x.className=cls;if(text!==undefined)x.textContent=text;return x;};
let experiment=null,guide=null,indices=[],active='static',generation=0;
const preloads=new Map(),openTables=new Set();
const sendHeight=()=>{if(window.parent!==window)window.parent.postMessage({type:'foundations-figure-height',id:config.id,height:Math.ceil(document.body.getBoundingClientRect().height)},location.origin);};
new ResizeObserver(sendHeight).observe(document.body);

// Keep feedback visible when switching to the approved fallback after a failure.
const toolbar=document.querySelector('.toolbar'),status=$('status');toolbar.after(status);
const actions=el('div','figure-actions');actions.append($('reset'));
const permalink=el('a','permalink','Open figure ↗');permalink.target='_blank';permalink.rel='noopener';permalink.setAttribute('aria-label','Open this figure with the current parameters');actions.append(permalink);toolbar.append(actions);
const chapterURL=new URL(`../../${config.id.split('--')[0]}/`,location.href);
if(config.bookLabel)chapterURL.hash=config.bookLabel.replace(/['’]/g,'-').replace(/-+/g,'-');
if(window.parent===window){
  document.body.classList.add('standalone');
  const header=el('header','figure-heading'),back=el('a','chapter-link','← Read in the chapter');back.href=chapterURL;
  header.append(back,el('h1',null,config.title));$('experiment').prepend(header);
}
const intro=el('p','exploration-prompt');$('controls').before(intro);
const presets=el('div','presets');presets.setAttribute('aria-label','Suggested comparisons');$('controls').after(presets);
const observation=el('aside','observation');observation.append(el('strong',null,'What to look for'),el('p'));$('metrics').after(observation);
document.querySelector('.grid-note').hidden=true;
const audio=el('audio');audio.controls=true;audio.preload='none';audio.src='../data/bird.wav';audio.setAttribute('aria-label','Play the recorded birdsong');audio.hidden=true;$('panels').after(audio);
$('static').tabIndex=0;

function updateLink(){
  const query=new URLSearchParams();
  if(experiment)experiment.controls.forEach((c,i)=>query.set(c.id,c.values[indices[i]]));
  if(active==='static'&&config.data)query.set('view','book');
  const target=location.pathname+(query.size?'?'+query:'');
  history.replaceState(null,'',target);permalink.href=new URL(target,location.href);
}
function choose(view){
  active=view;generation++;
  if(view==='static'){$('explore').setAttribute('aria-busy','false');document.body.dataset.rendering='false';if(!status.classList.contains('error'))status.textContent='';audio.pause();}
  $('explore').hidden=view!=='explore';$('static').hidden=view!=='static';
  for(const [id,v] of [['explore-tab','explore'],['static-tab','static']]){
    $(id).setAttribute('aria-selected',String(view===v));$(id).tabIndex=view===v?0:-1;
  }
  $('reset').hidden=view!=='explore';
  if(view==='explore'&&experiment)redraw();
  if(experiment||!config.data)updateLink();sendHeight();
}
$('explore-tab').onclick=()=>choose('explore');$('static-tab').onclick=()=>choose('static');
for(const tab of [$('explore-tab'),$('static-tab')])tab.onkeydown=e=>{
  if(!['ArrowLeft','ArrowRight','Home','End'].includes(e.key))return;
  e.preventDefault();const next=e.key==='Home'?$('explore-tab'):e.key==='End'?$('static-tab'):tab.id==='explore-tab'?$('static-tab'):$('explore-tab');
  if(!next.hidden){next.click();next.focus();}
};
$('reset').onclick=()=>{indices=experiment.controls.map(c=>c.default);syncControls();redraw();};
function syncControls(){
  experiment.controls.forEach((c,j)=>{
    const input=$('parameter-'+c.id),value=displayValue(c.values[indices[j]]);input.value=indices[j];
    if(c.type==='range'){$('value-'+c.id).textContent=value;input.setAttribute('aria-valuetext',value);}
  });
  for(const button of presets.querySelectorAll('button')){
    const values=guide.presets[Number(button.dataset.preset)].values;
    button.setAttribute('aria-pressed',String(Object.entries(values).every(([id,v])=>experiment.controls.some((c,j)=>c.id===id&&c.values[indices[j]]===v))));
  }
  updateLink();
}
function makeControls(){
  const saved=new URLSearchParams(location.search);
  indices=experiment.controls.map(c=>{const index=c.values.findIndex(x=>String(x)===saved.get(c.id));return index>=0?index:c.default;});
  $('controls').replaceChildren();
  experiment.controls.forEach((c,j)=>{
    const wrap=el('div','control'),label=el('label',null,c.label);label.htmlFor='parameter-'+c.id;let input;
    if(c.type==='range'){
      const output=el('output');output.id='value-'+c.id;output.htmlFor='parameter-'+c.id;label.append(output);
      input=el('input');input.type='range';input.min=0;input.max=c.values.length-1;input.step=1;
    }else{
      input=el('select');c.values.forEach((v,i)=>{const option=el('option',null,displayValue(v));option.value=i;input.append(option);});
    }
    input.id='parameter-'+c.id;
    input.addEventListener('input',()=>{indices[j]=Number(input.value);syncControls();redraw();});wrap.append(label,input);
    if(c.type==='range'){
      const marks=el('div','range-labels');marks.setAttribute('aria-hidden','true');
      for(const i of [0,Math.floor((c.values.length-1)/2),c.values.length-1]){const mark=el('span',null,displayValue(c.values[i]));mark.style.left=`${100*i/(c.values.length-1)}%`;marks.append(mark);}
      wrap.append(marks);
    }
    const explanation=controlHelp(experiment,c);
    if(explanation){const help=el('p','control-help',explanation);help.id='help-'+c.id;input.setAttribute('aria-describedby',help.id);wrap.append(help);}
    $('controls').append(wrap);
  });
  guide=getGuide(experiment);intro.textContent=guide.prompt;observation.querySelector('p').textContent=guide.observe;
  presets.replaceChildren(el('span',null,'Try'));
  guide.presets.forEach((preset,i)=>{
    const button=el('button',null,preset.label);button.dataset.preset=i;
    button.onclick=()=>{for(const [id,value] of Object.entries(preset.values)){const j=experiment.controls.findIndex(c=>c.id===id);indices[j]=experiment.controls[j].values.indexOf(value);}syncControls();redraw();};
    presets.append(button);
  });
  syncControls();
}
function tablePanel(panel){
  const wrap=el('div','table-scroll'),table=el('table'),head=el('thead'),row=el('tr');
  panel.columns.forEach(v=>{const th=el('th',null,v);th.scope='col';row.append(th);});head.append(row);table.append(head);
  const body=el('tbody');panel.rows.forEach(values=>{const tr=el('tr');values.forEach(v=>tr.append(el('td',null,v)));body.append(tr);});table.append(body);wrap.append(table);return wrap;
}
function plotPanel(panel,panelIndex){
  const wrap=el('div'),canvas=el('canvas');canvas.setAttribute('role','img');canvas.tabIndex=0;
  canvas.setAttribute('aria-label',`${panel.title}. Horizontal axis: ${panel.xlabel||'x'}. Vertical axis: ${panel.ylabel||'y'}. Arrow keys inspect values: left and right select a point, up and down select a series.`);wrap.append(canvas);
  const legend=el('div','legend');
  panel.series.forEach((s,i)=>{
    if(!s.label)return;
    const item=el('span'),mark=el('i',s.points?'point':s.dash?'dashed':'');mark.style.setProperty('--color',plotColor(s.color||colors[i%colors.length],document.documentElement.classList.contains('dark')));
    item.append(mark,document.createTextNode(s.label));legend.append(item);
  });
  if(legend.children.length>8){const more=el('details','legend-details');more.append(el('summary',null,`Series legend (${legend.children.length})`),legend);wrap.append(more);}else wrap.append(legend);
  const detail=el('div','plot-tooltip','Point to a curve, or focus the plot and use the arrow keys.');wrap.append(detail);
  const announce=el('span','sr-only');announce.setAttribute('aria-live','polite');wrap.append(announce);
  requestAnimationFrame(()=>{
    if(!canvas.isConnected)return;
    const transforms=drawPlot(canvas,panel,Math.max(240,wrap.clientWidth));
    const backdrop=document.createElement('canvas');backdrop.width=canvas.width;backdrop.height=canvas.height;backdrop.getContext('2d').drawImage(canvas,0,0);
    const ctx=canvas.getContext('2d');let selectedSeries=0,selectedPoint=0;
    const clear=()=>{ctx.save();ctx.setTransform(1,0,0,1,0,0);ctx.clearRect(0,0,canvas.width,canvas.height);ctx.drawImage(backdrop,0,0);ctx.restore();};
    function inspect(si,pi,keyboard=false){
      const series=panel.series[si];if(!series?.x.length)return;
      selectedSeries=si;selectedPoint=Math.max(0,Math.min(pi,series.x.length-1));
      const i=selectedPoint,px=transforms.X(series.x[i]),py=transforms.Y(series.y[i]),box=transforms.box;
      const text=`${series.label||`Series ${si+1}`}: ${panel.xlabel||'x'} = ${format(series.x[i])}, ${panel.ylabel||'y'} = ${format(series.y[i])}`;
      detail.textContent=text;if(keyboard)announce.textContent=text;clear();
      if(px<box.x||px>box.x+box.w||py<box.y||py>box.y+box.h||!Number.isFinite(px+py))return;
      ctx.save();ctx.strokeStyle=plotColor(series.color||colors[si%colors.length],document.documentElement.classList.contains('dark'));ctx.lineWidth=1;ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(px,box.y);ctx.lineTo(px,box.y+box.h);ctx.stroke();ctx.setLineDash([]);ctx.fillStyle=ctx.strokeStyle;ctx.beginPath();ctx.arc(px,py,4,0,2*Math.PI);ctx.fill();ctx.restore();
    }
    canvas.onpointermove=e=>{
      const r=canvas.getBoundingClientRect(),px=(e.clientX-r.left)/r.width*transforms.width,py=(e.clientY-r.top)/r.height*transforms.height;let best=null;
      panel.series.forEach((s,si)=>{for(let i=0;i<s.x.length;i++){const d=Math.hypot(transforms.X(s.x[i])-px,transforms.Y(s.y[i])-py);if(Number.isFinite(d)&&(!best||d<best.d))best={d,si,i};}});
      if(best)inspect(best.si,best.i);
    };
    canvas.onpointerleave=()=>{if(document.activeElement!==canvas){clear();detail.textContent='Point to a curve, or focus the plot and use the arrow keys.';}};
    canvas.onfocus=()=>inspect(selectedSeries,selectedPoint,true);
    canvas.onkeydown=e=>{
      if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','Home','End'].includes(e.key))return;
      e.preventDefault();let si=selectedSeries,pi=selectedPoint;
      if(e.key==='ArrowLeft')pi--;if(e.key==='ArrowRight')pi++;if(e.key==='Home')pi=0;if(e.key==='End')pi=panel.series[si].x.length-1;
      if(e.key==='ArrowUp')si=(si-1+panel.series.length)%panel.series.length;if(e.key==='ArrowDown')si=(si+1)%panel.series.length;
      inspect(si,pi,true);
    };
  });
  const data=el('details','chart-table');data.dataset.panel=panelIndex;data.append(el('summary',null,'Numerical values'));
  const fill=()=>{if(data.querySelector('table'))return;const rows=[];panel.series.forEach(s=>s.x.forEach((x,i)=>rows.push([s.label,format(x),format(s.y[i])])));data.append(tablePanel({columns:['Series',panel.xlabel||'x',panel.ylabel||'y'],rows}));};
  data.addEventListener('toggle',()=>{if(!data.isConnected)return;if(data.open){openTables.add(panelIndex);fill();}else openTables.delete(panelIndex);});
  if(openTables.has(panelIndex)){data.open=true;fill();}wrap.append(data);return wrap;
}

function imagePanel(panel){const wrap=el('div','image-wrap');wrap.style.setProperty('--ratio',`${panel.width}/${panel.height}`);const img=el('img');img.src=new URL('../data/'+panel.src,location.href);img.alt=panel.title;img.width=panel.width;img.height=panel.height;img.onload=sendHeight;wrap.append(img);
  const marks=[...(panel.boundaries||[]),...(panel.boxes||[])];if(marks.length||(panel.vectors||[]).length){const svg=document.createElementNS('http://www.w3.org/2000/svg','svg');svg.setAttribute('viewBox',`0 0 ${panel.width} ${panel.height}`);marks.forEach(b=>{const rect=document.createElementNS(svg.namespaceURI,'rect');for(const [k,v]of Object.entries({x:b.x,y:b.y,width:b.w||.1,height:b.h||.1,fill:'none',stroke:b.color,'stroke-width':Math.max(1,panel.width/180)}))rect.setAttribute(k,v);svg.append(rect);});for(const v of panel.vectors||[]){const a=document.createElementNS(svg.namespaceURI,'path');const x=v.x+v.dx,y=v.y+v.dy,angle=Math.atan2(v.dy,v.dx),r=.18;const d=`M ${v.x} ${v.y} L ${x} ${y} M ${x-r*Math.cos(angle-.6)} ${y-r*Math.sin(angle-.6)} L ${x} ${y} L ${x-r*Math.cos(angle+.6)} ${y-r*Math.sin(angle+.6)}`;a.setAttribute('d',d);a.setAttribute('fill','none');a.setAttribute('stroke','#e62421');a.setAttribute('stroke-width','.065');svg.append(a);}wrap.append(svg);}const outer=el('div');outer.append(wrap);if(panel.extent){const [x0,x1,y1,y0]=panel.extent;outer.append(el('p','image-coordinates',`${panel.xlabel||'x'}: ${format(x0)} to ${format(x1)} · ${panel.ylabel||'y'}: ${format(y0)} to ${format(y1)}`));}if(panel.probability){const key=el('div','probability-key');key.append(el('span',null,'0'),el('i'),el('span',null,'1'));outer.append(key,el('p','image-coordinates','Probability of class 1'));}if(panel.colorLabels){const leg=el('div','legend');for(const label of panel.colorLabels){const item=el('span'),mark=el('i');mark.style.setProperty('--color',label.color);item.append(mark,document.createTextNode(label.label));leg.append(item);}outer.append(leg);}return outer;}
function meshPanel(panel){const wrap=el('div','image-wrap');wrap.style.setProperty('--ratio','1');const canvas=el('canvas');canvas.width=600;canvas.height=600;canvas.setAttribute('role','img');canvas.setAttribute('aria-label',panel.title+' with triangles following the curved boundary');wrap.append(canvas);const c=canvas.getContext('2d');const bg=new Image();bg.src=new URL('../data/'+panel.background,location.href);bg.onload=()=>{const [x0,x1,y1,y0]=panel.extent||[0,1,1,0],X=x=>(x-x0)/(x1-x0)*600,Y=y=>(y-y0)/(y1-y0)*600;c.globalAlpha=.4;c.drawImage(bg,x0*bg.width,y0*bg.height,(x1-x0)*bg.width,(y1-y0)*bg.height,0,0,600,600);c.globalAlpha=1;c.strokeStyle='#146c78';c.lineWidth=.65;c.beginPath();for(const line of panel.lines){line.forEach(([x,y],i)=>i?c.lineTo(X(x),Y(y)):c.moveTo(X(x),Y(y)));}c.stroke();c.strokeStyle='#b94b44';c.lineWidth=2;c.beginPath();c.ellipse(X(.5),Y(.5),.28/(x1-x0)*600,.34/(y1-y0)*600,0,0,2*Math.PI);c.stroke();};return wrap;}


function preloadImage(src){
  const url=new URL('../data/'+src,location.href).href;
  if(!preloads.has(url)){
    const promise=new Promise((resolve,reject)=>{const image=new Image();image.onload=()=>resolve();image.onerror=()=>reject(new Error('Missing experiment image'));image.src=url;});
    preloads.set(url,promise);promise.catch(()=>preloads.delete(url));
  }
  return preloads.get(url);
}
async function redraw(){
  if(!experiment||active!=='explore')return;
  const ticket=++generation,key=indices.join(','),state=experiment.states[key];
  if(!state){fail(new Error('This parameter combination is unavailable.'));return;}
  $('explore').setAttribute('aria-busy','true');document.body.dataset.rendering='true';
  try{
    const images=state.panels.flatMap(p=>p.kind==='image'?[p.src]:p.kind==='mesh'?[p.background]:[]);
    if(images.some(src=>!preloads.has(new URL('../data/'+src,location.href).href)))status.textContent='Updating figure…';
    await Promise.all(images.map(preloadImage));
    if(ticket!==generation||active!=='explore')return;
    const frag=document.createDocumentFragment();
    state.panels.forEach((p,i)=>{const fig=el('figure','panel');fig.append(el('h3',null,p.title),p.kind==='plot'?plotPanel(p,i):p.kind==='table'?tablePanel(p):p.kind==='mesh'?meshPanel(p):imagePanel(p));frag.append(fig);});
    $('panels').replaceChildren(frag);
    $('metrics').replaceChildren(...(state.metrics||[]).map(m=>{const x=el('div','metric',m.label);x.append(el('strong',null,m.value));return x;}));
    $('metrics').hidden=!(state.metrics||[]).length;$('note').textContent=state.note||'';
    audio.hidden=!(config.id.startsWith('shannon--')&&state.parameters.source==='Birdsong');if(audio.hidden)audio.pause();
    status.textContent='';status.classList.remove('error');delete document.body.dataset.error;
    document.body.dataset.state=key;updateLink();sendHeight();
  }catch(error){if(ticket===generation)fail(error);}
  finally{if(ticket===generation){$('explore').setAttribute('aria-busy','false');document.body.dataset.rendering='false';}}
}
function fail(error){
  choose('static');status.classList.add('error');status.replaceChildren(document.createTextNode('The interactive view could not load. The book figure is shown below. '));
  const retry=el('button','retry','Try again');retry.onclick=()=>{experiment?choose('explore'):loadExperiment();};status.append(retry);
  document.body.dataset.error=String(error);document.body.dataset.rendering='false';sendHeight();
}
async function loadExperiment(){
  status.classList.remove('error');status.textContent='Preparing experiment…';
  const savedView=new URLSearchParams(location.search).get('view');
  try{
    const response=await fetch(config.data);if(!response.ok)throw new Error(`Experiment data: ${response.status}`);
    experiment=await response.json();makeControls();choose(savedView==='book'?'static':'explore');
    if(savedView==='book')status.textContent='';document.body.dataset.ready='interactive';delete document.body.dataset.error;
  }catch(error){fail(error);}
}
function setTheme(dark){if(document.documentElement.classList.contains('dark')===dark)return;document.documentElement.classList.toggle('dark',dark);if(experiment&&active==='explore')redraw();}
if(window.parent!==window){try{setTheme(window.parent.document.documentElement.classList.contains('dark'));}catch{}}
else setTheme(matchMedia('(prefers-color-scheme: dark)').matches);
window.addEventListener('message',e=>{if(e.origin===location.origin&&e.source===window.parent&&e.data?.type==='foundations-theme')setTheme(Boolean(e.data.dark));});
if(!config.data){$('explore-tab').hidden=true;$('reset').hidden=true;status.textContent='';choose('static');document.body.dataset.ready='static';}
else {choose('explore');loadExperiment();}
let previousWidth=window.innerWidth;
window.addEventListener('resize',()=>{if(window.innerWidth===previousWidth)return;previousWidth=window.innerWidth;if(experiment&&active==='explore')redraw();});
