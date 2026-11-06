"use strict";
const test=require("node:test");
const assert=require("node:assert/strict");
const {transcriptProgress,filterTranscripts,transcriptPage,polymeraseMotion,transcriptColor}=require("../static/transcription-view.js");

test("zero-speed polymerases are distinct from biochemical pauses and released RNA",()=>{
  const moving={state:"nascent",paused:false,elongation_rate:0.5};
  const stalled={state:"nascent",paused:false,elongation_rate:0};
  const paused={state:"nascent",paused:true,elongation_rate:0};
  assert.equal(polymeraseMotion(moving),"elongating");
  assert.equal(polymeraseMotion(stalled),"stalled");
  assert.equal(polymeraseMotion(paused),"paused");
  assert.equal(polymeraseMotion({state:"mature",paused:false,elongation_rate:0}),"released");
  assert.notEqual(transcriptColor(stalled),transcriptColor(paused));
  assert.notEqual(transcriptColor(stalled),transcriptColor(moving));
});

test("RNA length includes both gene endpoint bins and preserves fractional growth",()=>{
  const gene={start:10,end:14};
  assert.equal(transcriptProgress({length_nt:0},gene,200),0);
  assert.equal(transcriptProgress({length_nt:250.5},gene,200),0.2505);
  assert.equal(transcriptProgress({length_nt:1000},gene,200),1);
  assert.equal(transcriptProgress({length_nt:1200},gene,200),1);
  assert.equal(transcriptProgress({length_nt:-1},gene,200),0);
});

test("transcript search keeps gene and RNA-state boundaries",()=>{
  const records=[
    {id:"rna-1",polymerase_id:"pol-1",gene_id:"a",state:"nascent"},
    {id:"rna-2",polymerase_id:"pol-2",gene_id:"a",state:"mature"},
    {id:"rna-3",polymerase_id:"pol-3",gene_id:"b",state:"mature"}
  ];
  assert.deepEqual(filterTranscripts(records,"a","mature").map(t=>t.id),["rna-2"]);
  assert.deepEqual(filterTranscripts(records,"a",""," POL-2 ").map(t=>t.id),["rna-2"]);
  assert.deepEqual(filterTranscripts(records,"a","nascent","rna-2"),[]);
  assert.deepEqual(filterTranscripts(records,"a","","rna-3"),[]);
  assert.equal(records.length,3);
});

test("all 5000 allowed molecules remain reachable through bounded table pages",()=>{
  const records=Array.from({length:5000},(_,id)=>({id:`rna-${id}`}));
  const seen=new Set();
  for(let page=0;page<200;page++){
    const result=transcriptPage(records,page);
    assert.equal(result.rows.length,25);
    result.rows.forEach(record=>seen.add(record.id));
  }
  assert.equal(seen.size,5000);
  assert.equal(transcriptPage(records,999).page,199);
  assert.equal(transcriptPage(records,-1).page,0);
  assert.deepEqual(transcriptPage([],200),{page:0,pages:1,total:0,start:0,rows:[]});
});
