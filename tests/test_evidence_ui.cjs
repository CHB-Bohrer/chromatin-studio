"use strict";
// Pure helpers of the Evidence view (static/evidence-view.js). Run with: node --test tests/
const test=require("node:test"),assert=require("node:assert/strict");
const E=require("../static/evidence-view.js");
const esc=value=>String(value).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const s=(label,x,y,extra={})=>({label,role:"model",style:"points",group:null,x,y,...extra});

test("numbers are compact, zero stays a value and missing values are a dash",()=>{
  assert.equal(E.fmt(0),"0");assert.equal(E.fmt(0.0297),"0.0297");assert.equal(E.fmt(1.5e-5),"1.5e-5");
  assert.equal(E.fmt(123456),"1.23e5");assert.equal(E.fmt(1e-12),"1e-12");
  for(const value of [null,undefined,NaN,Infinity,"3"])assert.equal(E.fmt(value),"—");
});
test("a log axis ignores zero and underflow instead of stretching to 1e-300",()=>{
  const [low,high]=E.domain([1,0.5,0,1e-300,0.001],{log:true});
  assert.ok(low>1e-4&&low<0.001,String(low));assert.ok(high>1&&high<1.5,String(high));   // padding: 5% of the decades
  assert.deepEqual(E.domain([],{log:true}),[1,10]);
  const [zeroLow]=E.domain([0.2,3],{zero:true});assert.equal(zeroLow,0);
  assert.ok(E.logTicks(0.001,1).includes(0.01));
});
test("groups share a hue, measurements outside a group stay neutral, slots keep a fixed order",()=>{
  const list=[s("Engine, A",[1],[1],{group:"A"}),s("Measured, A",[1],[1],{group:"A",role:"experiment"}),s("Engine, B",[1],[2],{group:"B"}),s("Published",[1],[3],{role:"experiment"}),s("Closed form",[1],[4],{role:"reference",style:"line"})];
  const colors=E.colorPlan(list);
  assert.equal(colors[0],E.SLOTS[0]);assert.equal(colors[1],E.SLOTS[0]);assert.equal(colors[2],E.SLOTS[1]);
  assert.notEqual(colors[3],E.SLOTS[2]);assert.notEqual(colors[4],E.SLOTS[2]);
  assert.deepEqual(E.slotPlan(list),[0,0,1,-1,-1]);
});
test("from the fourth coloured group on, engine markers and legend glyphs become squares",()=>{
  const figure={id:"f",title:"Four",kind:"xy",x_label:"x",y_label:"y",x_scale:"linear",y_scale:"linear",bands:[],
    series:[0,1,2,3].map(k=>s(`Engine ${k}`,[1,2],[k+1,k+2]))};
  const model=E.drawFigure(figure,{width:400,height:260,escape:esc});
  assert.equal((model.svg.match(/<rect x=/g)||[]).length-1,2,"two squares for the fourth series (plus the clip rectangle)");
  const legend=E.legendItems(figure);
  assert.match(legend[3].html,/<rect/);assert.doesNotMatch(legend[0].html,/<rect/);
  assert.equal(E.SAFE_HUES,3);
});
test("points below the log axis are hidden and counted, labels are escaped",()=>{
  const figure={id:"f",title:"Survival <tail>",kind:"xy",x_label:"t",y_label:"P",x_scale:"linear",y_scale:"log",bands:[],
    series:[s("Engine",[0,1,2,3],[1,0.1,0,1e-300],{style:"line"})]};
  const model=E.drawFigure(figure,{width:400,height:260,escape:esc});
  assert.equal(model.hidden,2);assert.equal(model.points.length,2);
  assert.match(model.svg,/Survival &lt;tail&gt;/);assert.doesNotMatch(model.svg,/<tail>/);
  assert.match(model.svg,/role="img"/);
});
test("category figures keep every series of a column together for the tooltip",()=>{
  const figure={id:"c",title:"Ratios",kind:"category",categories:["Reb1","Pho4"],x_label:"factor",y_label:"ratio",x_scale:"linear",y_scale:"log",
    bands:[{label:"Measured",role:"experiment",category:"Pho4",x:null,y:[40,null]}],
    series:[s("Engine",["Reb1","Pho4"],[0.03,0.01]),s("Measured",["Reb1","Pho4"],[0.02,0.009],{role:"experiment"})]};
  const model=E.drawFigure(figure,{width:400,height:260,escape:esc});
  assert.equal(model.points.filter(p=>p.category==="Reb1").length,2);
  assert.ok(model.layout.category);
  const legend=E.legendItems(figure);assert.equal(legend.length,3);
});
test("tables show a range only when its ends differ after formatting",()=>{
  const figure={id:"t",title:"T",kind:"xy",x_label:"x",y_label:"y",series:[s("Engine",[1,2],[0.0297,0.5],{low:[0.0297,0.4],high:[0.0297,0.6]})]};
  const table=E.tableRows(figure);
  assert.equal(table.rows[0][3],"—");assert.equal(table.rows[1][3],"0.4 to 0.6");
});
test("check strips use the stated window, an open end and a log scale for wide spans",()=>{
  const open=E.checkGeometry({rule:"inside stated interval",low:40,high:null,observed:120,expected:40});
  assert.deepEqual(open.accept,[40,null]);assert.equal(open.log,false);
  const wide=E.checkGeometry({rule:"relative tolerance",expected:1,allowed:0.1,observed:53,sem:2});
  assert.equal(wide.log,true);
  assert.equal(E.expectedText({rule:"inside stated interval",low:0.5,high:null}),"≥ 0.5");
  assert.equal(E.expectedText({rule:"relative tolerance",expected:2,allowed:0.1}),"2 ± 0.1");
});
test("the comparison basis is read only from the start of a note",()=>{
  assert.equal(E.basis("Tuned: the scanned hazard closest to the measured ratio is selected."),"tuned");
  assert.equal(E.basis("Prediction of rebinding among methylated neighbours"),"prediction");
  assert.equal(E.basis("Input-matched: the release factor was set from this ratio"),"input-matched");
  assert.equal(E.basis("Uses the tuned hazard"),null);assert.equal(E.basis(null),null);
});
test("the overview prefers the latest standard report and sums checks by evidence level",()=>{
  const reports=[{study:"a",scale:"smoke",analyzed_utc:"2026-11-21T10:00"},{study:"a",scale:"standard",analyzed_utc:"2026-11-15T10:00"},
    {study:"a",scale:"standard",analyzed_utc:"2026-11-20T10:00"},{study:"b",error:"unreadable"}];
  const best=E.latestByStudy(reports);
  assert.equal(best.get("a").analyzed_utc,"2026-11-20T10:00");assert.equal(best.has("b"),false);
  const sum=E.totals([{check_summary:{total:3,passed:2,by_evidence:{numerical_check:{total:2,passed:2},experimental_reference:{total:1,passed:0}}}}]);
  assert.equal(sum.total,3);assert.equal(sum.passed,2);assert.equal(sum.levels.experimental_reference.total,1);
});
