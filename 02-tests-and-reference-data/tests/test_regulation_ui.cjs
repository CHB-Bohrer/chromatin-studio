"use strict";
const test=require("node:test");
const assert=require("node:assert/strict");
const {uniqueMolecules,regionInterval,localMarkSummary,localNucleosomeSummary}=require("../static/regulation-view.js");

test("regulatory circuit totals count shared TFs and coactivators once",()=>{
  const regions=[{tf_molecules:[1,2],recruited_molecules:[4,5]},{tf_molecules:[2,3],recruited_molecules:[5]},{}];
  assert.equal(uniqueMolecules(regions,"tf_molecules"),3);
  assert.equal(uniqueMolecules(regions,"recruited_molecules"),2);
  assert.equal(uniqueMolecules([],"recruited_molecules"),0);
});

test("one coarse nucleosome unit remains a full 200 bp interval",()=>{
  assert.deepEqual(regionInterval({start:3,end:3},200),{start:600,end:800});
  assert.deepEqual(regionInterval({start:398,end:399},200),{start:79600,end:80000});
  assert.deepEqual(regionInterval({start:30,end:31},10000),{start:300000,end:320000});
});

test("overlapping elements do not double-count local marked beads",()=>{
  const config={marks:[{id:"active",name:"Active",color:"#abc"},{id:"silent",name:"Silent",color:"#def"}]};
  const snapshot={marks:{active:[false,true,true,false,false],silent:[false,false,true,true,false]}};
  const result=localMarkSummary(config,snapshot,[{start:1,end:2},{start:2,end:3}],null);
  assert.equal(result.beads,3);
  assert.equal(result.marked,3);
  assert.deepEqual(result.marks.map(m=>m.count),[2,2]);
  assert.equal(snapshot.marks.active[1],true);
});

test("genes without elements report marks over their own complete interval",()=>{
  const result=localMarkSummary({marks:[{id:"m"}]},{marks:{m:[true,false,true]}},[],{start:1,end:2});
  assert.equal(result.beads,2);assert.equal(result.marked,1);assert.equal(result.marks[0].count,1);
  assert.deepEqual(localMarkSummary({marks:[]},{},[],null),{beads:0,marked:0,marks:[]});
});

test("overlapping elements count each actual occupied or vacant bead once",()=>{
  const snapshot={nucleosomes:{enabled:true,occupied:[1,0,1,0]}};
  assert.deepEqual(localNucleosomeSummary(snapshot,[{start:0,end:2},{start:1,end:3}],null),
    {enabled:true,beads:4,occupied:2,vacant:2});
  assert.deepEqual(localNucleosomeSummary(snapshot,[],{start:1,end:1}),
    {enabled:true,beads:1,occupied:0,vacant:1});
});

test("disabled or unavailable occupancy is never reported as vacant",()=>{
  assert.deepEqual(localNucleosomeSummary({nucleosomes:{enabled:false,occupied:[]}},[],{start:0,end:2}),
    {enabled:false,beads:0,occupied:null,vacant:null});
  assert.deepEqual(localNucleosomeSummary({nucleosomes:{enabled:true,occupied:[]}},[],{start:0,end:2}),
    {enabled:true,beads:0,occupied:null,vacant:null});
});

test("causal map excludes disabled elements, deduplicates shared molecules, and separates promoter from nascent state",()=>{
  const {causalSummary}=require("../static/regulation-view.js");
  const config={genes:[{id:"g",start:0,end:5}],marks:[{id:"m"}],regulatory_regions:[{id:"p",kind:"promoter",target_gene:"g",start:0,end:1},{id:"e",kind:"enhancer",target_gene:"g",start:1,end:2},{id:"disabled",kind:"enhancer",target_gene:"g",enabled:false,start:3,end:4}]};
  const snapshot={genes:[{id:"g",on:false,nascent:2}],marks:{m:[false,true,false,true]},nucleosomes:{enabled:true,occupied:[1,0,1,0,0]},regulation:{genes:[{id:"g",mode:"regions"}],regions:[{id:"p",tf_molecules:[1],recruited_molecules:[3],contact:true},{id:"e",tf_molecules:[1,2],recruited_molecules:[3],contact:false},{id:"disabled",tf_molecules:[99],recruited_molecules:[88],contact:true}]}};
  const summary=causalSummary(config,snapshot,"g");
  assert.equal(summary.tf,2);assert.equal(summary.coactivators,1);assert.equal(summary.contacts,0);assert.equal(summary.totalEnhancers,1);
  assert.equal(summary.access.beads,3);assert.equal(summary.access.vacant,1);assert.equal(summary.marks.marked,1);
  assert.equal(summary.promoterOn,false);assert.equal(summary.nascentPositive,true);
});
test("unavailable rate reports are not displayed as measured zero occupancy",()=>{
  const {causalSummary}=require("../static/regulation-view.js"),summary=causalSummary({genes:[{id:"g",start:0,end:1}]},{},"g");
  assert.equal(summary.tf,null);assert.equal(summary.coactivators,null);assert.equal(summary.contacts,null);assert.equal(summary.promoterOn,null);assert.equal(summary.nascentPositive,null);
});
