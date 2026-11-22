/* Evidence view: saved component-study reports, drawn from the figure data each
   study returns. Read-only: it never runs a simulation, applies a draft or
   touches the live run. A study changes its plots in its own analyze(); this
   file only draws them, with one generic renderer. */
"use strict";
(function(global){
  const finite=value=>typeof value==="number"&&Number.isFinite(value);
  const LEVELS=[
    {id:"numerical_check",title:"Numerical checks",text:"The engine reproduces its own stated rule."},
    {id:"analytical_reference",title:"Analytical references",text:"The engine reproduces a closed-form result."},
    {id:"experimental_reference",title:"Experimental comparisons",text:"An engine output is compared with a published measurement under the stated mapping. Inside the range is not calibration or validation."}];
  // Engine groups take these hues in fixed order (validated on the dark chart surface).
  // Text never wears them; identity also comes from the legend and the mark shape.
  const SLOTS=["#3987e5","#d95926","#199e70","#c98500","#d55181","#0ca30c","#9085e9","#e66767"];
  const NEUTRAL={experiment:"#e8edf3",reference:"#a5b2c1",model:"#a5b2c1"};
  const SURFACE="#111923";
  const STATUS={pass:"#0ca30c",fail:"#d03b3b"};
  const round=value=>Math.round(value*10)/10;

  function fmt(value,digits=3){
    if(!finite(value))return "—";
    if(value===0)return "0";
    const size=Math.abs(value);
    if(size>=1e5||size<1e-3)return value.toExponential(digits-1).replace(/\.?0+e/,"e").replace("e+","e");
    return String(Number(value.toPrecision(digits)));
  }
  function niceStep(span,count){const raw=span/Math.max(1,count),power=10**Math.floor(Math.log10(raw)),f=raw/power;return (f<1.5?1:f<3?2:f<7?5:10)*power;}
  function linearTicks(low,high,count=5){
    if(!finite(low)||!finite(high))return [];
    if(high<=low)return [low];
    const step=niceStep(high-low,count),ticks=[];
    for(let value=Math.ceil(low/step-1e-9)*step;value<=high+step*1e-9;value+=step)ticks.push(Math.abs(value)<step*1e-9?0:Number(value.toPrecision(12)));
    return ticks;
  }
  function logTicks(low,high){
    if(!(low>0)||!(high>low))return [];
    const first=Math.floor(Math.log10(low)),last=Math.ceil(Math.log10(high)),decades=last-first,every=Math.max(1,Math.ceil(decades/7)),ticks=[];
    for(let power=first;power<=last;power++){
      if((power-first)%every)continue;
      for(const m of decades<=2?[1,2,5]:[1]){const value=Number((m*10**power).toPrecision(12));if(value>=low*(1-1e-9)&&value<=high*(1+1e-9))ticks.push(value);}
    }
    return ticks.length>=2?ticks:linearTicks(low,high,3).filter(value=>value>0);
  }
  // Data domain with padding; a log axis ignores values <= 0. "zero" starts a
  // non-negative linear axis at 0 when the data come close to it.
  function domain(values,{log=false,zero=false}={}){
    let list=values.filter(value=>finite(value)&&(!log||value>0));
    if(!list.length)return log?[1,10]:[0,1];
    if(log){const top=Math.max(...list);list=list.filter(value=>value>=top*1e-12);}   // underflow is not data
    let low=Math.min(...list),high=Math.max(...list);
    if(log){if(low===high){low/=2;high*=2;}const pad=Math.max(1.12,(high/low)**.05);return [low/pad,high*pad];}
    if(zero&&low>=0&&low<=.35*high)low=0;
    if(low===high){const d=Math.abs(low)||1;return [low-d*.5,high+d*.5];}
    const pad=(high-low)*.06;return [low===0?0:low-pad,high+pad];
  }
  function scaler([r0,r1],[d0,d1],log){
    const a=log?Math.log10(d0):d0,b=log?Math.log10(d1):d1;
    return value=>{if(!finite(value)||(log&&value<=0))return null;return r0+((log?Math.log10(value):value)-a)/(b-a)*(r1-r0);};
  }
  // Series in one group share a hue; only groups with an engine series take a hue.
  // Ungrouped measurements are near-white, ungrouped references gray.
  function slotPlan(seriesList){
    const slot=new Map(),key=s=>s.group!=null?"g:"+s.group:"s:"+s.label;
    for(const s of seriesList){if(s.role!=="model")continue;const k=key(s);if(!slot.has(k)&&slot.size<SLOTS.length)slot.set(k,slot.size);}
    return seriesList.map(s=>slot.has(key(s))?slot.get(key(s)):-1);
  }
  function colorPlan(seriesList){return slotPlan(seriesList).map((k,i)=>k>=0?SLOTS[k]:NEUTRAL[seriesList[i].role]||NEUTRAL.reference);}
  // Only the first three hues stay apart for every pair under colour-vision
  // deficiency (validated against the chart surface); later ones only from their
  // neighbours. Studies keep to three coloured groups per figure; if one does not,
  // engine markers from the fourth group on become squares, so identity never
  // rests on hue alone.
  const SAFE_HUES=3;
  function squareSlot(k){return k>=SAFE_HUES;}
  function wrap(text,width=16,lines=2){
    const out=[];let line="";
    for(const word of String(text).split(/\s+/)){if(!line)line=word;else if((line+" "+word).length<=width)line+=" "+word;else{out.push(line);line=word;}}
    if(line)out.push(line);
    if(out.length>lines){out.length=lines;out[lines-1]=out[lines-1].slice(0,Math.max(1,width-1))+"…";}
    return out;
  }
  function usableSeries(figure){return (figure?.series||[]).filter(s=>s&&Array.isArray(s.x)&&Array.isArray(s.y)&&s.x.length===s.y.length);}
  function linePath(coords,step){
    let d="",open=false,last=null;
    for(const c of coords){
      if(!c){open=false;continue;}
      if(!open){d+=`M${round(c[0])} ${round(c[1])}`;open=true;}
      else d+=step?`H${round(c[0])}V${round(c[1])}`:`L${round(c[0])} ${round(c[1])}`;
      last=c;
    }
    return last?d:"";
  }
  function marker(role,[x,y],color,square){
    if(role==="experiment")return `<path d="M${round(x)} ${round(y-6)}L${round(x+6)} ${round(y)}L${round(x)} ${round(y+6)}L${round(x-6)} ${round(y)}Z" fill="${SURFACE}" stroke="${color}" stroke-width="2"/>`;
    if(role==="reference")return `<circle cx="${round(x)}" cy="${round(y)}" r="3.5" fill="${SURFACE}" stroke="${color}" stroke-width="1.5"/>`;
    if(square)return `<rect x="${round(x-5)}" y="${round(y-5)}" width="10" height="10" rx="1.5" fill="${color}" stroke="${SURFACE}" stroke-width="2"/>`;
    return `<circle cx="${round(x)}" cy="${round(y)}" r="5" fill="${color}" stroke="${SURFACE}" stroke-width="2"/>`;
  }
  let clipCounter=0;
  // Pure: figure data + pixel size -> SVG markup and the hover points (in SVG units).
  function drawFigure(figure,{width=640,height=320,escape=value=>String(value)}={}){
    const category=figure.kind==="category",categories=category?(figure.categories||[]).map(String):[];
    const margin={left:66,right:20,top:16,bottom:category?60:50};
    const x0=margin.left,x1=Math.max(x0+60,width-margin.right),y1=margin.top,y0=Math.max(y1+60,height-margin.bottom);
    const xlog=!category&&figure.x_scale==="log",ylog=figure.y_scale==="log";
    const list=usableSeries(figure),colors=colorPlan(list),slots=slotPlan(list),bands=(figure.bands||[]).filter(b=>b&&(Array.isArray(b.x)||Array.isArray(b.y)));
    const yValues=[],xValues=[];
    list.forEach(s=>s.y.forEach((value,i)=>{if(!finite(value))return;yValues.push(value,s.low?.[i],s.high?.[i]);if(!category)xValues.push(s.x[i]);}));
    bands.forEach(b=>{(b.y||[]).forEach(value=>yValues.push(value));if(!category)(b.x||[]).forEach(value=>xValues.push(value));});
    const yDomain=domain(yValues,{log:ylog,zero:!ylog}),xDomain=category?null:domain(xValues,{log:xlog,zero:!xlog});
    const Y=scaler([y0,y1],yDomain,ylog),column=category?(x1-x0)/Math.max(1,categories.length):0;
    const categoryX=name=>{const i=categories.indexOf(String(name));return i<0?null:x0+column*(i+.5);};
    const X=category?categoryX:scaler([x0,x1],xDomain,xlog);
    const active=list.map((s,i)=>s.y.some(finite)?i:-1).filter(i=>i>=0);
    const dodge=category?Math.min(16,column*.7/Math.max(1,active.length)):0;
    const offset=i=>category?(active.indexOf(i)-(active.length-1)/2)*dodge:0;
    const id="ev-clip-"+(++clipCounter),axis=[],grid=[];
    (ylog?logTicks(...yDomain):linearTicks(...yDomain,5)).forEach(value=>{const y=Y(value);if(y==null||y<y1-.5||y>y0+.5)return;grid.push(`<line class="ev-grid" x1="${x0}" x2="${x1}" y1="${round(y)}" y2="${round(y)}"/>`);axis.push(`<text class="ev-tick" x="${x0-8}" y="${round(y)+4}" text-anchor="end">${escape(fmt(value))}</text>`);});
    if(category)categories.forEach((name,i)=>{const cx=x0+column*(i+.5);wrap(name,Math.max(8,Math.floor(column/6.6))).forEach((line,k)=>axis.push(`<text class="ev-tick" x="${round(cx)}" y="${y0+18+k*13}" text-anchor="middle">${escape(line)}</text>`));});
    else (xlog?logTicks(...xDomain):linearTicks(...xDomain,Math.max(3,Math.round((x1-x0)/110)))).forEach(value=>{const x=X(value);if(x==null||x<x0-.5||x>x1+.5)return;axis.push(`<line class="ev-axis" x1="${round(x)}" x2="${round(x)}" y1="${y0}" y2="${y0+5}"/><text class="ev-tick" x="${round(x)}" y="${y0+19}" text-anchor="middle">${escape(fmt(value))}</text>`);});
    axis.push(`<line class="ev-axis" x1="${x0}" x2="${x1}" y1="${y0}" y2="${y0}"/>`);
    axis.push(`<text class="ev-axis-title" x="${round((x0+x1)/2)}" y="${height-8}" text-anchor="middle">${escape(figure.x_label||"")}${xlog?" · log scale":""}</text>`);
    axis.push(`<text class="ev-axis-title" transform="translate(15 ${round((y0+y1)/2)}) rotate(-90)" text-anchor="middle">${escape(figure.y_label||"")}${ylog?" · log":""}</text>`);
    const clamp=(value,lo,hi)=>Math.min(hi,Math.max(lo,value)),shaded=[];
    bands.forEach(b=>{
      const color=b.role==="experiment"?NEUTRAL.experiment:NEUTRAL.reference;
      let left=x0,right=x1,top=y1,bottom=y0;
      if(Array.isArray(b.y)){const [lo,hi]=b.y;top=finite(hi)?(Y(hi)??y1):y1;bottom=finite(lo)?(Y(lo)??y0):y0;}
      if(category&&b.category!=null){const cx=categoryX(b.category);if(cx==null)return;left=cx-column*.42;right=cx+column*.42;}
      else if(!category&&Array.isArray(b.x)){const [lo,hi]=b.x;left=finite(lo)?(X(lo)??x0):x0;right=finite(hi)?(X(hi)??x1):x1;}
      top=clamp(top,y1,y0);bottom=clamp(bottom,y1,y0);left=clamp(left,x0,x1);right=clamp(right,x0,x1);
      if(bottom-top<.5||right-left<.5)return;
      shaded.push(`<rect x="${round(left)}" y="${round(top)}" width="${round(right-left)}" height="${round(bottom-top)}" fill="${color}" fill-opacity=".1"/>`);
      if(Array.isArray(b.y))[[b.y[1],top],[b.y[0],bottom]].forEach(([edge,y])=>{if(finite(edge))shaded.push(`<line x1="${round(left)}" x2="${round(right)}" y1="${round(y)}" y2="${round(y)}" stroke="${color}" stroke-opacity=".4"/>`);});
    });
    const lines={experiment:[],reference:[],model:[]},marks={model:[],experiment:[],reference:[]},errors=[],points=[];let hidden=0;
    list.forEach((s,i)=>{
      const color=colors[i],role=lines[s.role]?s.role:"model",dx=offset(i);
      const inside=(value,dom,log)=>!log||value>=dom[0]*(1-1e-9);
      const coords=s.x.map((xv,k)=>{const px=X(xv),py=Y(s.y[k]);return px==null||py==null||!inside(s.y[k],yDomain,ylog)||(!category&&!inside(xv,xDomain,xlog))?null:[px+dx,py];});
      hidden+=s.y.filter((value,k)=>finite(value)&&!coords[k]).length;
      if(!category&&s.style!=="points"){const d=linePath(coords,s.style==="step");if(d)lines[role].push(`<path d="${d}" class="ev-line-${role}" stroke="${color}"/>`);}
      const marked=category||s.style==="points"||s.style==="connected";
      coords.forEach((c,k)=>{
        if(!c)return;
        const lo=s.low?.[k],hi=s.high?.[k];
        if(finite(lo)&&finite(hi)&&hi>lo){const a=Y(hi)??y1,b=Y(lo)??y0;errors.push(`<path d="M${round(c[0])} ${round(a)}V${round(b)}M${round(c[0]-4)} ${round(a)}H${round(c[0]+4)}M${round(c[0]-4)} ${round(b)}H${round(c[0]+4)}" class="ev-error" stroke="${color}"/>`);}
        if(marked)marks[role].push(marker(role,c,color,squareSlot(slots[i])));
        points.push({x:c[0],y:c[1],series:i,label:s.label,role,color,xv:s.x[k],yv:s.y[k],lo:finite(lo)?lo:null,hi:finite(hi)?hi:null,category:category?String(s.x[k]):null});
      });
    });
    const label=`${figure.title||"Figure"}. ${list.length} series${bands.length?` and ${bands.length} shaded ranges`:""}. The values are listed in the table below the chart.`;
    const svg=`<svg class="ev-svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img" aria-label="${escape(label)}"><defs><clipPath id="${id}"><rect x="${x0-8}" y="${y1-8}" width="${x1-x0+16}" height="${y0-y1+16}"/></clipPath></defs>${grid.join("")}${axis.join("")}<g clip-path="url(#${id})">${shaded.join("")}${lines.experiment.join("")}${lines.reference.join("")}${lines.model.join("")}${errors.join("")}${marks.model.join("")}${marks.experiment.join("")}${marks.reference.join("")}</g><circle class="ev-focus" r="9" cx="0" cy="0" visibility="hidden"/></svg>`;
    return {svg,points,hidden,colors,layout:{width,height,x0,x1,y0,y1,column,category,categories}};
  }
  function glyph(kind,color,square){
    const svg=body=>`<svg class="ev-glyph" viewBox="0 0 24 12" aria-hidden="true">${body}</svg>`;
    if(square&&kind==="model")return svg(`<rect x="7.5" y="1.5" width="9" height="9" rx="1.5" fill="${color}" stroke="${SURFACE}" stroke-width="1.5"/>`);
    if(square&&kind==="model-connected")return svg(`<path d="M2 6H22" stroke="${color}" stroke-width="2"/><rect x="8" y="2" width="8" height="8" rx="1.5" fill="${color}" stroke="${SURFACE}" stroke-width="1.5"/>`);
    if(kind==="band-experiment"||kind==="band-reference")return svg(`<rect x="2" y="1" width="20" height="10" rx="2" fill="${color}" fill-opacity=".18" stroke="${color}" stroke-opacity=".5"/>`);
    if(kind==="experiment-line")return svg(`<path d="M2 6H22" stroke="${color}" stroke-width="5" stroke-opacity=".35" stroke-linecap="round"/>`);
    if(kind==="experiment")return svg(`<path d="M12 1L17 6L12 11L7 6Z" fill="${SURFACE}" stroke="${color}" stroke-width="2"/>`);
    if(kind==="reference-line")return svg(`<path d="M2 6H22" stroke="${color}" stroke-width="2" stroke-dasharray="5 3"/>`);
    if(kind==="reference")return svg(`<circle cx="12" cy="6" r="3.5" fill="${SURFACE}" stroke="${color}" stroke-width="1.5"/>`);
    if(kind==="model-line")return svg(`<path d="M2 6H22" stroke="${color}" stroke-width="2"/>`);
    if(kind==="model-connected")return svg(`<path d="M2 6H22" stroke="${color}" stroke-width="2"/><circle cx="12" cy="6" r="4" fill="${color}" stroke="${SURFACE}" stroke-width="1.5"/>`);
    return svg(`<circle cx="12" cy="6" r="4.5" fill="${color}" stroke="${SURFACE}" stroke-width="1.5"/>`);
  }
  function glyphKind(s,category){
    const line=!category&&s.style!=="points"&&s.style!=="connected";
    if(s.role==="experiment")return line?"experiment-line":"experiment";
    if(s.role==="reference")return line?"reference-line":"reference";
    return line?"model-line":s.style==="connected"&&!category?"model-connected":"model";
  }
  function legendItems(figure){
    const list=usableSeries(figure),colors=colorPlan(list),slots=slotPlan(list),items=[],seen=new Set(),category=figure.kind==="category";
    list.forEach((s,i)=>{const key="s"+s.label;if(seen.has(key))return;seen.add(key);items.push({label:s.label,html:glyph(glyphKind(s,category),colors[i],s.role==="model"&&squareSlot(slots[i]))});});
    (figure.bands||[]).forEach(b=>{const key="b"+b.label+b.role;if(seen.has(key))return;seen.add(key);const color=b.role==="experiment"?NEUTRAL.experiment:NEUTRAL.reference;items.push({label:b.label,html:glyph("band-"+(b.role==="experiment"?"experiment":"reference"),color)});});
    return items;
  }
  function tableRows(figure){
    const list=usableSeries(figure);
    if(figure.kind==="category"){
      const categories=(figure.categories||[]).map(String);
      return {head:[figure.x_label||"category",...list.map(s=>s.label)],rows:categories.map(name=>[name,...list.map(s=>{const k=s.x.map(String).indexOf(name);if(k<0||!finite(s.y[k]))return "—";const lo=s.low?.[k],hi=s.high?.[k];return fmt(s.y[k])+(finite(lo)&&finite(hi)&&fmt(hi)!==fmt(lo)?` (${fmt(lo)} to ${fmt(hi)})`:"");})])};
    }
    const rows=[];
    list.forEach(s=>s.x.forEach((x,k)=>{if(!finite(s.y[k]))return;const lo=s.low?.[k],hi=s.high?.[k];rows.push([s.label,fmt(x),fmt(s.y[k]),finite(lo)&&finite(hi)&&fmt(hi)!==fmt(lo)?`${fmt(lo)} to ${fmt(hi)}`:"—"]);}));
    return {head:["series",figure.x_label||"x",figure.y_label||"y","±1 SEM range"],rows};
  }
  // Check strips: acceptance window (shaded), expected value (tick), observation ± 1 SEM.
  function checkGeometry(check){
    const interval=check.rule==="inside stated interval";
    const accept=interval?[finite(check.low)?check.low:null,finite(check.high)?check.high:null]
      :finite(check.expected)&&finite(check.allowed)?[check.expected-check.allowed,check.expected+check.allowed]:[null,null];
    const values=[check.observed,check.expected,...accept];
    if(finite(check.observed)&&finite(check.sem))values.push(check.observed-check.sem,check.observed+check.sem);
    const known=values.filter(finite);
    const log=known.length>1&&known.every(value=>value>0)&&Math.max(...known)/Math.min(...known)>=20;
    return {interval,accept,log,domain:domain(values,{log})};
  }
  function checkStrip(check,{width=240,height=34}={}){
    const g=checkGeometry(check),x0=8,x1=width-8,X=scaler([x0,x1],g.domain,g.log),mid=height/2,parts=[`<line class="ev-strip-track" x1="${x0}" x2="${x1}" y1="${mid}" y2="${mid}"/>`];
    const [lo,hi]=g.accept;
    if(finite(lo)||finite(hi)){const a=finite(lo)?(X(lo)??x0):x0,b=finite(hi)?(X(hi)??x1):x1,left=Math.max(x0,Math.min(a,b)),right=Math.min(x1,Math.max(a,b));parts.push(`<rect class="ev-strip-accept" x="${round(left)}" y="${mid-8}" width="${round(Math.max(2,right-left))}" height="16" rx="3"/>`);}
    if(finite(check.expected)){const e=X(check.expected);if(e!=null)parts.push(`<line class="ev-strip-expected" x1="${round(e)}" x2="${round(e)}" y1="${mid-10}" y2="${mid+10}"/>`);}
    if(finite(check.observed)){const o=X(check.observed);if(o!=null){
      if(finite(check.sem)&&check.sem>0){const a=X(check.observed-check.sem),b=X(check.observed+check.sem);if(a!=null&&b!=null)parts.push(`<line class="ev-strip-sem" x1="${round(a)}" x2="${round(b)}" y1="${mid}" y2="${mid}"/>`);}
      parts.push(`<circle cx="${round(o)}" cy="${mid}" r="5" fill="${check.passed?STATUS.pass:STATUS.fail}" stroke="${SURFACE}" stroke-width="2"/>`);}}
    if(g.log)parts.push(`<text class="ev-strip-scale" x="${x1}" y="${height-1}" text-anchor="end">log</text>`);
    return `<svg class="ev-strip" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" aria-hidden="true">${parts.join("")}</svg>`;
  }
  function expectedText(check){
    if(check.rule==="inside stated interval"){
      const lo=check.low,hi=check.high;
      if(finite(lo)&&finite(hi))return `${fmt(lo)} to ${fmt(hi)}`;
      if(finite(lo))return `≥ ${fmt(lo)}`;
      if(finite(hi))return `≤ ${fmt(hi)}`;
    }
    if(finite(check.expected)&&finite(check.allowed))return `${fmt(check.expected)} ± ${fmt(check.allowed)}`;
    return fmt(check.expected);
  }
  // How an experimental comparison was set up, when its note says so.
  function basis(note){
    const match=/^\s*(tuned|prediction|predicted|input-matched|input multiplier|input)\b/i.exec(String(note||""));
    if(!match)return null;
    const word=match[1].toLowerCase();
    return word.startsWith("tuned")?"tuned":word.startsWith("predict")?"prediction":"input-matched";
  }
  // One report per study for the overview: the latest standard run, else the latest of any scale.
  function latestByStudy(reports){
    const best=new Map(),rank=r=>[(r.scale==="standard")?1:0,r.analyzed_utc||""];
    for(const r of reports||[]){
      if(!r||r.error||!r.study)continue;
      const current=best.get(r.study);
      if(!current){best.set(r.study,r);continue;}
      const [a,b]=[rank(r),rank(current)];
      if(a[0]>b[0]||(a[0]===b[0]&&a[1]>b[1]))best.set(r.study,r);
    }
    return best;
  }
  function totals(reports){
    const levels={};LEVELS.forEach(level=>levels[level.id]={total:0,passed:0});
    let total=0,passed=0;
    for(const r of reports){const s=r.check_summary||{};total+=s.total||0;passed+=s.passed||0;for(const level of LEVELS){const b=s.by_evidence?.[level.id];if(b){levels[level.id].total+=b.total||0;levels[level.id].passed+=b.passed||0;}}}
    return {total,passed,levels};
  }
  const helpers={finite,fmt,linearTicks,logTicks,domain,slotPlan,colorPlan,drawFigure,legendItems,tableRows,checkGeometry,checkStrip,expectedText,basis,latestByStudy,totals,SLOTS,SAFE_HUES};
  if(typeof module!=="undefined"&&module.exports)module.exports=helpers;
  if(typeof document==="undefined")return;

  const esc=value=>String(value??"").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const day=value=>typeof value==="string"&&value.length>=10?value.slice(0,10):"—";
  class EvidenceView{
    constructor(container,{session={},onNavigate=()=>{}}={}){
      this.root=container;this.session=session;this.onNavigate=onNavigate;this.active=true;
      session.reports??={};
      this.click=e=>{const target=e.target.closest("[data-ev-action],[data-ev-select]");if(!target||!this.root.contains(target))return;if(target.dataset.evSelect){this.select(target.dataset.evSelect);return;}const action=target.dataset.evAction;if(action==="back")this.onNavigate("system");if(action==="reload"){this.session.index=null;this.session.reports={};this.load();}};
      this.change=e=>{if(e.target.dataset.evVersion!==undefined)this.select(e.target.value,{keepStudy:true});};
      container.addEventListener("click",this.click);container.addEventListener("change",this.change);
      this.resize=typeof ResizeObserver==="function"?new ResizeObserver(()=>{cancelAnimationFrame(this.frame);this.frame=requestAnimationFrame(()=>this.drawFigures());}):null;
      this.render();if(!session.index)this.load();
    }
    async request(path){const response=await fetch(path);const data=await response.json();if(!response.ok)throw new Error(data.error||"The evidence request failed.");return data;}
    async load(){
      this.session.loading=true;this.session.error="";this.render();
      try{
        this.session.index=await this.request("/api/evidence");
        const best=latestByStudy(this.session.index.reports);
        const keep=this.session.index.reports.some(r=>r.path===this.session.selected);
        if(!keep){const first=(this.session.index.catalog||[]).map(study=>best.get(study.id)).find(Boolean)||[...best.values()][0];this.session.selected=first?.path||null;}
      }catch(error){this.session.error=error.message;}
      finally{this.session.loading=false;}
      if(this.session.selected)await this.select(this.session.selected,{keepStudy:true});else this.render();
    }
    async select(path,{keepStudy=false}={}){
      if(!path)return;
      if(!keepStudy){const entry=this.session.index?.reports?.find(r=>r.path===path);if(entry)this.session.study=entry.study;}
      this.session.selected=path;
      if(!this.session.reports[path]){
        this.session.busy=path;this.render();
        try{this.session.reports[path]=await this.request("/api/evidence/report?path="+encodeURIComponent(path));}
        catch(error){this.session.error=error.message;}
        finally{this.session.busy=null;}
      }
      if(this.active){this.render();this.root.querySelector("#ev-report-heading")?.focus({preventScroll:true});}
    }
    tiles(best){
      const reports=[...best.values()],sum=totals(reports),catalog=this.session.index?.catalog||[];
      const tile=(label,value,detail)=>`<div class="ev-tile"><span>${esc(label)}</span><strong>${value}</strong><small>${esc(detail)}</small></div>`;
      return `<section class="ev-tiles" aria-label="Summary of the latest report per study">${tile("Studies with a report",`${reports.length}<em> of ${catalog.length}</em>`,"latest standard run per study")}${tile("Checks inside their rule",`${sum.passed}<em> of ${sum.total}</em>`,"all evidence levels")}${LEVELS.map(level=>tile(level.title,`${sum.levels[level.id].passed}<em> of ${sum.levels[level.id].total}</em>`,level.id==="experimental_reference"?"published measurements":level.id==="analytical_reference"?"closed-form results":"engine rules")).join("")}</section>`;
    }
    listHTML(best){
      const index=this.session.index,catalog=index?.catalog||[],byStudy=new Map();
      (index?.reports||[]).forEach(r=>{if(!r.study)return;if(!byStudy.has(r.study))byStudy.set(r.study,[]);byStudy.get(r.study).push(r);});
      const current=this.session.reports[this.session.selected]?.study||this.session.study;
      const row=(study,report)=>{
        if(!report)return `<div class="ev-row ev-row-empty"><span class="ev-row-title">${esc(study.title)}</span><span class="ev-row-meta">Not run on this machine yet</span></div>`;
        const s=report.check_summary||{},total=s.total||0,passed=s.passed||0,flags=[];
        if(report.missing_jobs||report.partial_results)flags.push(`partial: ${report.complete_jobs??"?"} of ${report.expected_jobs??"?"} runs`);
        if(report.changed_since_run?.length)flags.push("code changed since run");
        if(!report.figures)flags.push("no figure data");
        return `<button class="ev-row" data-ev-select="${esc(report.path)}" aria-current="${study.id===current?"true":"false"}"><span class="ev-row-title">${esc(study.title||report.title)}</span><span class="ev-row-meta">${esc(report.scale||"")} · ${esc(day(report.analyzed_utc))} · ${passed} of ${total} checks</span><span class="ev-meter" aria-hidden="true"><i class="ev-meter-pass" style="width:${total?passed/total*100:0}%"></i><i class="ev-meter-fail" style="width:${total?(total-passed)/total*100:0}%"></i></span>${flags.length?`<span class="ev-row-flags">${esc(flags.join(" · "))}</span>`:""}</button>`;
      };
      const known=new Set(catalog.map(study=>study.id));
      const extra=[...best.values()].filter(r=>!known.has(r.study)).map(r=>row({id:r.study,title:r.title||r.study},r));
      const broken=(index?.reports||[]).filter(r=>r.error);
      return `<nav class="ev-list" aria-label="Component studies">${catalog.map(study=>row(study,best.get(study.id))).join("")}${extra.length?`<h4>Other reports</h4>${extra.join("")}`:""}${broken.length?`<p class="ev-small">${broken.length} report file${broken.length>1?"s":""} could not be read: ${broken.map(r=>esc(r.path)).join(", ")}.</p>`:""}</nav>`;
    }
    detailHTML(){
      const path=this.session.selected,report=this.session.reports[path];
      if(this.session.busy)return `<p class="ev-status" role="status">Loading the report…</p>`;
      if(!report)return `<div class="ev-empty"><h3>Select a study</h3><p>Each report lists every check that ran, with its figures and sources.</p></div>`;
      const entry=this.session.index?.reports?.find(r=>r.path===path)||{},versions=(this.session.index?.reports||[]).filter(r=>r.study===report.study&&!r.error);
      const folder=p=>p.split("/").slice(0,-2).join("/")||"top level";
      const history=[...versions.map(r=>({path:r.path,label:`${r.scale||""} run · folder ${folder(r.path)} · analyzed ${day(r.analyzed_utc)}`})),...versions.flatMap(r=>(r.history||[]).map(h=>({path:h,label:`earlier analysis · folder ${folder(h.replace(/\/report-history\/[^/]+$/,"/x"))} · ${h.split("/").pop().replace(/^report-|\.json$/g,"")}`})))];
      const chips=[];
      chips.push(report.source_matches_current?`<span class="ev-chip ev-chip-ok"><b aria-hidden="true">✓</b>Simulated with the current code</span>`:`<span class="ev-chip ev-chip-warn"><b aria-hidden="true">!</b>Code changed since the run${report.changed_since_run?.length?": "+esc(report.changed_since_run.join(", ")):""}</span>`);
      if(finite(report.expected_jobs))chips.push(`<span class="ev-chip ${report.complete_jobs<report.expected_jobs?"ev-chip-warn":""}">${report.complete_jobs??"?"} of ${report.expected_jobs} runs complete${report.partial_results?` · ${report.partial_results} partial`:""}</span>`);
      else chips.push(`<span class="ev-chip">${report.results_used??"?"} results used</span>`);
      chips.push(`<span class="ev-chip">Analyzed ${esc(day(report.analyzed_utc))}</span>`);
      let html=`<article class="ev-report"><header class="ev-report-head"><span class="ev-kicker">${esc(report.study)} · ${esc(report.scale||"")} scale</span><h3 id="ev-report-heading" tabindex="-1">${esc(report.title||report.study)}</h3>${report.question?`<p class="ev-question">${esc(report.question)}</p>`:""}<div class="ev-chips">${chips.join("")}</div>${history.length>1?`<label class="ev-version">Report<select data-ev-version>${history.map(h=>`<option value="${esc(h.path)}" ${h.path===path?"selected":""}>${esc(h.label)}</option>`).join("")}</select></label>`:""}</header>`;
      if((report.findings||[]).length)html+=`<section class="ev-findings" aria-label="Findings"><h4>What the runs show</h4><ul>${report.findings.map(f=>`<li>${esc(f)}</li>`).join("")}</ul></section>`;
      const figures=report.figures||[];
      if(figures.length)html+=`<section class="ev-figures" aria-label="Figures"><div class="ev-section-head"><h4>Figures</h4><div class="ev-role-key" aria-label="How marks are encoded"><span>${glyph("model","#b8c7d6")}Engine output</span><span>${glyph("experiment","#e8edf3")}Published measurement</span><span>${glyph("reference-line","#a5b2c1")}Closed form or model input</span><span>${glyph("band-experiment","#e8edf3")}Measured or accepted range</span></div></div><div class="ev-figure-grid">${figures.map((figure,i)=>this.figureHTML(figure,i)).join("")}</div></section>`;
      else html+=`<p class="ev-notice">${entry.format<2||report.format<2?"This report was written before studies returned figure data. ":""}Re-analyze its folder with the current code to draw its plots: <code>python -m chromatin_studio.components analyze &lt;study folder&gt;</code>. The checks below are unaffected.</p>`;
      html+=this.checksHTML(report.checks||[]);
      if((report.sources||[]).length)html+=`<section class="ev-sources" aria-label="Sources"><h4>Sources and reported values</h4><ul>${report.sources.map(s=>`<li><a href="${esc(/^https:\/\//.test(s.url||"")?s.url:"#")}" target="_blank" rel="noopener noreferrer">${esc(s.citation)}</a><span>${esc(s.measured)}</span></li>`).join("")}</ul></section>`;
      const scope=Array.isArray(report.scope)?report.scope:report.scope?[report.scope]:[];
      if(scope.length)html+=`<details class="ev-scope"><summary>Scope and interpretation</summary><ul>${scope.map(s=>`<li>${esc(typeof s==="string"?s:JSON.stringify(s))}</li>`).join("")}</ul></details>`;
      return html+`</article>`;
    }
    figureHTML(figure,i){
      const items=legendItems(figure),table=tableRows(figure);
      return `<figure class="ev-figure"><figcaption><h5>${esc(figure.title||figure.id)}</h5></figcaption><div class="ev-plot" data-ev-plot="${i}" tabindex="0" role="group" aria-label="${esc((figure.title||figure.id)+". Arrow keys step through the values.")}" aria-describedby="ev-note-${i}"><div class="ev-plot-slot"></div><div class="ev-tooltip" aria-live="polite" hidden></div></div>${items.length>1?`<div class="ev-legend">${items.map(item=>`<span>${item.html}${esc(item.label)}</span>`).join("")}</div>`:""}<p class="ev-note" id="ev-note-${i}">${esc(figure.note||"")}<span class="ev-hidden-note" hidden></span></p><details class="ev-table"><summary>Values as a table</summary><div class="ev-table-wrap"><table><thead><tr>${table.head.map(h=>`<th>${esc(h)}</th>`).join("")}</tr></thead><tbody>${table.rows.map(row=>`<tr>${row.map((cell,k)=>k?`<td>${esc(cell)}</td>`:`<th>${esc(cell)}</th>`).join("")}</tr>`).join("")}</tbody></table></div></details></figure>`;
    }
    checksHTML(checks){
      const groups=LEVELS.map(level=>({level,checks:checks.filter(c=>c.evidence===level.id)})).filter(g=>g.checks.length);
      if(!groups.length)return `<p class="ev-notice">This report has no checks.</p>`;
      return `<section class="ev-checks" aria-label="Checks"><div class="ev-section-head"><h4>Every check that ran</h4><p class="ev-small">Shaded: accepted window. Tick: expected value. Dot: observation (green pass, red fail; the label says which), with ±1 SEM.</p></div>${groups.map(({level,checks:list})=>`<div class="ev-check-group"><h5>${esc(level.title)} <span>${list.filter(c=>c.passed).length} of ${list.length} inside</span></h5><p class="ev-small">${esc(level.text)}</p>${list.map(c=>{const tag=c.evidence==="experimental_reference"?basis(c.note):null;return `<div class="ev-check"><span class="ev-verdict ${c.passed?"pass":"fail"}"><b aria-hidden="true">${c.passed?"✓":"✕"}</b>${c.passed?"Inside":"Outside"}</span><div class="ev-check-text"><strong>${esc(c.label)}</strong>${tag?`<em class="ev-tag">${esc(tag)}</em>`:""}<span>${esc(fmt(c.observed))}${finite(c.sem)?` ± ${esc(fmt(c.sem))}`:""} vs ${esc(expectedText(c))}${c.unit&&c.unit!=="dimensionless"?" "+esc(c.unit):""}</span>${c.note?`<small>${esc(c.note)}</small>`:""}</div>${checkStrip(c)}</div>`;}).join("")}</div>`).join("")}</section>`;
    }
    render(){
      if(!this.active)return;
      const s=this.session,index=s.index,best=latestByStudy(index?.reports);
      let html=`<div class="ev-page"><header class="ev-hero"><div><span class="ev-kicker">EVIDENCE · WHAT WAS ACTUALLY TESTED</span><h2>Component tests against references and experiments</h2><p>Each study isolates one mechanism and compares the engine with its own rule, a closed-form result or a published measurement. This page draws the saved reports; it runs nothing and leaves the live model untouched.</p></div><div class="ev-hero-actions"><button data-ev-action="reload" class="secondary-button">Reload reports</button><button data-ev-action="back" class="secondary-button">← Return to simulation</button></div></header><div class="ev-strip-notes"><span>Saved reports only</span><span>Three evidence levels kept apart</span><span>Reduced units unless a mapping is stated</span><span>Passing is not calibration or validation</span></div>`;
      if(s.error)html+=`<div class="ev-error" role="alert">${esc(s.error)}</div>`;
      if(!index){this.root.innerHTML=html+`<p class="ev-status" role="status">${s.loading?"Reading saved reports…":"The evidence index is unavailable."}</p></div>`;return;}
      if(!index.available)html+=`<section class="ev-connect"><h3>${index.configured?"The evidence folder was not found":"No evidence folder connected"}</h3><p>${esc(index.hint||"")}</p><pre>python -m chromatin_studio.server --evidence &lt;folder with study outputs&gt;</pre><p class="ev-small">The registered studies are listed below. Run one with <code>python -m chromatin_studio.components run &lt;study&gt; --output &lt;folder&gt;</code>.</p></section>`;
      html+=this.tiles(best)+`<div class="ev-layout">${this.listHTML(best)}<section class="ev-detail" aria-live="polite">${this.detailHTML()}</section></div></div>`;
      this.root.innerHTML=html;
      const detail=this.root.querySelector(".ev-detail");
      if(this.resize){this.resize.disconnect();if(detail)this.resize.observe(detail);}
      this.drawFigures();
    }
    drawFigures(){
      if(!this.active)return;
      const report=this.session.reports[this.session.selected];if(!report)return;
      this.root.querySelectorAll(".ev-plot").forEach(el=>{
        const figure=(report.figures||[])[Number(el.dataset.evPlot)],slot=el.querySelector(".ev-plot-slot");if(!figure||!slot)return;
        const width=Math.max(300,Math.min(980,Math.round(slot.clientWidth||640))),height=Math.round(Math.max(240,Math.min(340,width*.58)));
        if(el._ev&&el._ev.layout.width===width)return;
        const model=drawFigure(figure,{width,height,escape:esc});slot.innerHTML=model.svg;el._ev=model;
        const note=el.parentElement.querySelector(".ev-hidden-note");if(note){note.hidden=!model.hidden;note.textContent=model.hidden?` ${model.hidden} value${model.hidden>1?"s are":" is"} too small for the log axis (zero, or far below its range) and not drawn; the table lists every value.`:"";}
        if(!el._evBound){el._evBound=true;this.bindPlot(el);}
      });
    }
    bindPlot(el){
      const tip=el.querySelector(".ev-tooltip");let order=[],index=-1;
      const ring=()=>el.querySelector(".ev-focus");
      const hide=()=>{tip.hidden=true;ring()?.setAttribute("visibility","hidden");};
      const show=(anchor,rows,title)=>{
        const model=el._ev,svg=el.querySelector("svg");if(!model||!svg)return;
        const scale=svg.getBoundingClientRect().width/model.layout.width||1;
        tip.replaceChildren();
        if(title){const h=document.createElement("div");h.className="ev-tip-title";h.textContent=title;tip.append(h);}
        rows.forEach(p=>{const row=document.createElement("div");row.className="ev-tip-row";const key=document.createElement("i");key.style.background=p.color;const value=document.createElement("strong");value.textContent=fmt(p.yv)+(p.lo!=null&&p.hi!=null&&fmt(p.lo)!==fmt(p.hi)?` (${fmt(p.lo)} to ${fmt(p.hi)})`:"");const label=document.createElement("span");label.textContent=p.label+(model.layout.category?"":` · ${fmt(p.xv)}`);row.append(key,value,label);tip.append(row);});
        tip.hidden=false;
        const left=anchor.x*scale,top=anchor.y*scale;
        tip.style.left=Math.round(left)+"px";tip.style.top=Math.round(top)+"px";
        tip.classList.toggle("ev-tip-left",left>el.clientWidth*.6);
        const r=ring();if(r){r.setAttribute("cx",round(anchor.x));r.setAttribute("cy",round(anchor.y));r.setAttribute("visibility","visible");}
      };
      const nearest=(px,py)=>{
        const model=el._ev;if(!model)return null;const {layout,points}=model;
        if(layout.category){const i=Math.floor((px-layout.x0)/layout.column);const name=layout.categories[i];if(name==null||px<layout.x0||px>layout.x1)return null;const rows=points.filter(p=>p.category===name);if(!rows.length)return null;const anchor=rows.reduce((a,b)=>Math.abs(b.y-py)<Math.abs(a.y-py)?b:a);return {anchor,rows,title:name};}
        let best=null,distance=Infinity;for(const p of points){const d=(p.x-px)**2+(p.y-py)**2;if(d<distance){distance=d;best=p;}}
        return best&&distance<=34*34?{anchor:best,rows:[best],title:null}:null;
      };
      el.addEventListener("pointermove",e=>{const model=el._ev,svg=el.querySelector("svg");if(!model||!svg)return;const rect=svg.getBoundingClientRect(),k=model.layout.width/rect.width;const hit=nearest((e.clientX-rect.left)*k,(e.clientY-rect.top)*k);if(hit)show(hit.anchor,hit.rows,hit.title);else hide();});
      el.addEventListener("pointerleave",hide);
      el.addEventListener("focus",()=>{order=[...(el._ev?.points||[])].sort((a,b)=>a.x-b.x||a.y-b.y);index=order.length?0:-1;if(index>=0)show(order[0],[order[0]],order[0].category);});
      el.addEventListener("blur",hide);
      el.addEventListener("keydown",e=>{if(!order.length)return;const moves={ArrowRight:1,ArrowDown:1,ArrowLeft:-1,ArrowUp:-1};if(e.key==="Escape"){hide();return;}if(e.key==="Home")index=0;else if(e.key==="End")index=order.length-1;else if(moves[e.key])index=(index+moves[e.key]+order.length)%order.length;else return;e.preventDefault();const p=order[index];show(p,[p],p.category);});
    }
    destroy(){this.active=false;this.resize?.disconnect();cancelAnimationFrame(this.frame);this.root.removeEventListener("click",this.click);this.root.removeEventListener("change",this.change);this.root.innerHTML="";}
  }
  global.ChromatinEvidence={...helpers,render:(container,options)=>new EvidenceView(container,options)};
})(typeof window!=="undefined"?window:globalThis);
