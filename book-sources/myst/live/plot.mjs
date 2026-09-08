/** Small deterministic canvas renderer. Numeric axes always form a full box. */
export const colors=['#146c78','#b94b44','#4778a3','#c48a32','#78578e','#607b43','#ad597c'];
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
  const height=330,scale=window.devicePixelRatio||1;canvas.width=Math.round(width*scale);canvas.height=height*scale;canvas.style.aspectRatio=`${width}/${height}`;const c=canvas.getContext('2d');c.scale(scale,scale);c.fillStyle='white';c.fillRect(0,0,width,height);
  const box={x:62,y:12,w:width-78,h:height-63};let xb=panel.xlim||extent(panel.series.flatMap(s=>s.x),panel.logx),yb=panel.ylim||extent(panel.series.flatMap(s=>s.y),panel.logy);
  if(panel.series.some(s=>s.bars)&&!panel.ylim)yb[0]=Math.min(0,yb[0]);
  if(panel.equal){const unit=Math.max((xb[1]-xb[0])/box.w,(yb[1]-yb[0])/box.h);const mx=(xb[0]+xb[1])/2,my=(yb[0]+yb[1])/2;xb=[mx-unit*box.w/2,mx+unit*box.w/2];yb=[my-unit*box.h/2,my+unit*box.h/2];}
  const X=v=>coordinate(v,xb,box.x,box.w,panel.logx),Y=v=>coordinate(v,yb,box.y+box.h,-box.h,panel.logy);
  c.font='12px ui-sans-serif, sans-serif';c.lineWidth=.6;c.strokeStyle='#dae4e8';c.fillStyle='#506672';
  for(const x of ticks(xb,panel.logx)){const pos=X(x);c.beginPath();c.moveTo(pos,box.y);c.lineTo(pos,box.y+box.h);c.stroke();c.textAlign='center';c.fillText(format(x),pos,box.y+box.h+18);}
  for(const y of ticks(yb,panel.logy)){const pos=Y(y);c.beginPath();c.moveTo(box.x,pos);c.lineTo(box.x+box.w,pos);c.stroke();c.textAlign='right';c.fillText(format(y),box.x-7,pos+4);}
  c.strokeStyle='#728b97';c.lineWidth=1;c.strokeRect(box.x,box.y,box.w,box.h);c.font='14px ui-sans-serif, sans-serif';c.fillStyle='#203746';c.textAlign='center';c.fillText(panel.xlabel||'',box.x+box.w/2,height-6);c.save();c.translate(14,box.y+box.h/2);c.rotate(-Math.PI/2);c.fillText(panel.ylabel||'',0,0);c.restore();
  c.save();c.beginPath();c.rect(box.x,box.y,box.w,box.h);c.clip();
  panel.series.forEach((s,index)=>{const color=s.color||colors[index%colors.length];c.strokeStyle=color;c.fillStyle=color;c.lineWidth=s.width||1.8;c.setLineDash(s.dash?[5,4]:[]);let begun=false;c.beginPath();
    for(let i=0;i<s.x.length;i++){const x=s.x[i],y=s.y[i];if(!finite(x)||!finite(y)||(panel.logx&&x<=0)||(panel.logy&&y<=0)){begun=false;continue;}const px=X(x),py=Y(y);
      if(s.bars){const half=Math.max(1,Math.min(20,Math.abs(X((s.x[i+1]??x+1))-px)*.4));const zero=Y(panel.logy?yb[0]:0);c.globalAlpha=.65;c.fillRect(px-half,Math.min(py,zero),2*half,Math.abs(zero-py));c.globalAlpha=1;}
      else if(s.points){c.moveTo(px+(s.size||2.5),py);c.arc(px,py,s.size||2.5,0,2*Math.PI);}
      else if(s.stems){c.moveTo(px,Y(panel.logy?yb[0]:0));c.lineTo(px,py);}
      else if(!begun){c.moveTo(px,py);begun=true;}else c.lineTo(px,py);
    }if(s.points){c.globalAlpha=.7;c.fill();c.globalAlpha=1;}else c.stroke();
  });c.restore();
  return {X,Y,xb,yb,box,width,height};
}
