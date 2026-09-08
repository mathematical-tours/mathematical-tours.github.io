/** Small deterministic canvas renderer. Numeric axes always form a full box. */
export const colors=['#146c78','#b94b44','#4778a3','#c48a32','#78578e','#607b43','#ad597c'];
/** Keep series recognizable while lifting dark inks off the night background. */
export function plotColor(color,dark=false){
  if(!dark||!/^#[0-9a-f]{6}$/i.test(color))return color;
  return '#'+[1,3,5].map(i=>Math.round(parseInt(color.slice(i,i+2),16)*.64+255*.36).toString(16).padStart(2,'0')).join('');
}
export const finite=x=>typeof x==='number'&&Number.isFinite(x);
export function format(x){if(!finite(x))return '—';if(x===0)return '0';if(Math.abs(x)>=10000||Math.abs(x)<.001)return x.toExponential(1);return Number(x.toPrecision(3)).toString();}
export function extent(values,log=false){const a=values.filter(x=>finite(x)&&(!log||x>0));if(!a.length)return log?[.1,1]:[0,1];let lo=Math.min(...a),hi=Math.max(...a);if(lo===hi){const pad=Math.abs(lo)*.1||.5;lo-=pad;hi+=pad}if(log)return [Math.max(lo*.8,Number.MIN_VALUE),hi*1.2];const pad=(hi-lo)*.06;return [lo-pad,hi+pad];}
export function coordinate(value,bounds,start,length,log=false){const fn=log?Math.log10:x=>x;return start+length*(fn(value)-fn(bounds[0]))/(fn(bounds[1])-fn(bounds[0]));}
export function ticks(bounds,log=false){
  if(log){let out=[];for(let p=Math.ceil(Math.log10(bounds[0]));p<=Math.floor(Math.log10(bounds[1]));p++)out.push(10**p);return out.length?out:bounds;}
  const raw=(bounds[1]-bounds[0])/4,magnitude=10**Math.floor(Math.log10(raw));
  const step=[1,2,2.5,5,10].map(x=>x*magnitude).sort((a,b)=>Math.abs(Math.log(a/raw))-Math.abs(Math.log(b/raw)))[0];
  const start=Math.ceil(bounds[0]/step),end=Math.floor(bounds[1]/step);
  return Array.from({length:Math.max(0,end-start+1)},(_,i)=>Number(((start+i)*step).toPrecision(12)));
}
export function drawPlot(canvas,panel,width=470){
  const height=Math.min(350,Math.max(280,width*.72)),scale=window.devicePixelRatio||1;
  canvas.width=Math.round(width*scale);canvas.height=Math.round(height*scale);canvas.style.aspectRatio=`${width}/${height}`;
  const c=canvas.getContext('2d'),style=getComputedStyle(canvas),dark=document.documentElement.classList.contains('dark');
  const ink=style.getPropertyValue('--ink').trim()||'#203746',muted=style.getPropertyValue('--muted').trim()||'#536975',line=style.getPropertyValue('--line').trim()||'#cddce1';
  c.scale(scale,scale);c.fillStyle=style.getPropertyValue('--paper').trim()||'#fff';c.fillRect(0,0,width,height);
  const box={x:62,y:12,w:width-80,h:height-78};let xb=[...(panel.xlim||extent(panel.series.flatMap(s=>s.x),panel.logx))],yb=[...(panel.ylim||extent(panel.series.flatMap(s=>s.y),panel.logy))];
  if(panel.series.some(s=>s.bars)&&!panel.ylim&&!panel.logy)yb[0]=Math.min(0,yb[0]);
  if(panel.equal){const unit=Math.max((xb[1]-xb[0])/box.w,(yb[1]-yb[0])/box.h);const mx=(xb[0]+xb[1])/2,my=(yb[0]+yb[1])/2;xb=[mx-unit*box.w/2,mx+unit*box.w/2];yb=[my-unit*box.h/2,my+unit*box.h/2];}
  const X=v=>coordinate(v,xb,box.x,box.w,panel.logx),Y=v=>coordinate(v,yb,box.y+box.h,-box.h,panel.logy);
  c.font='12px ui-sans-serif, sans-serif';c.lineWidth=.6;c.strokeStyle=line;c.fillStyle=muted;
  const xt=ticks(xb,panel.logx),stride=Math.max(1,Math.ceil(xt.length/Math.max(2,box.w/58)));
  for(const [i,x] of xt.entries()){if(i%stride)continue;const pos=X(x);c.beginPath();c.moveTo(pos,box.y);c.lineTo(pos,box.y+box.h);c.stroke();c.textAlign='center';c.fillText(format(x),pos,box.y+box.h+18);}
  const yt=ticks(yb,panel.logy),ystride=Math.max(1,Math.ceil(yt.length/(box.h/32)));
  for(const [i,y] of yt.entries()){if(i%ystride)continue;const pos=Y(y);c.beginPath();c.moveTo(box.x,pos);c.lineTo(box.x+box.w,pos);c.stroke();c.textAlign='right';c.fillText(format(y),box.x-7,pos+4);}
  c.strokeStyle=muted;c.lineWidth=1;c.strokeRect(box.x,box.y,box.w,box.h);c.font='13px ui-sans-serif, sans-serif';c.fillStyle=ink;c.textAlign='center';
  // Wrap long axis labels at words instead of shrinking them into illegibility.
  function label(text,x,y,maxWidth){
    const lines=[''];for(const word of text.split(' ')){const i=lines.length-1,next=(lines[i]+' '+word).trim();if(lines[i]&&c.measureText(next).width>maxWidth)lines.push(word);else lines[i]=next;}
    lines.forEach((line,i)=>c.fillText(line,x,y+15*(i-(lines.length-1)/2)));
  }
  label(panel.xlabel||'',box.x+box.w/2,height-17,box.w+20);c.save();c.translate(14,box.y+box.h/2);c.rotate(-Math.PI/2);label(panel.ylabel||'',0,0,box.h);c.restore();
  c.save();c.beginPath();c.rect(box.x,box.y,box.w,box.h);c.clip();
  panel.series.forEach((s,index)=>{const color=plotColor(s.color||colors[index%colors.length],dark);c.strokeStyle=color;c.fillStyle=color;c.lineWidth=s.width||1.8;c.setLineDash(s.dash?[5,4]:[]);let begun=false;c.beginPath();
    for(let i=0;i<s.x.length;i++){const x=s.x[i],y=s.y[i];if(!finite(x)||!finite(y)||(panel.logx&&x<=0)||(panel.logy&&y<=0)){begun=false;continue;}const px=X(x),py=Y(y);
      if(s.bars){const spacing=s.barWidth??Math.abs((s.x[i+1]??x)-(s.x[i+1]===undefined?s.x[i-1]??x-1:x));const fraction=s.barWidth?0.5:0.4,left=X(x-spacing*fraction),right=X(x+spacing*fraction);const zero=Y(panel.logy?yb[0]:0);c.globalAlpha=.65;c.fillRect(left,Math.min(py,zero),Math.max(.5,right-left),Math.abs(zero-py));c.globalAlpha=1;}
      else if(s.points){c.moveTo(px+(s.size||2.5),py);c.arc(px,py,s.size||2.5,0,2*Math.PI);}
      else if(s.stems){c.moveTo(px,Y(panel.logy?yb[0]:0));c.lineTo(px,py);}
      else if(s.segments){if(i%2)c.lineTo(px,py);else c.moveTo(px,py);}
      else if(!begun){c.moveTo(px,py);begun=true;}else c.lineTo(px,py);
    }if(s.points){c.globalAlpha=.7;c.fill();c.globalAlpha=1;}else c.stroke();
  });c.restore();
  return {X,Y,xb,yb,box,width,height};
}
