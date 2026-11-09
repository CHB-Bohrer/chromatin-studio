"use strict";
const test=require("node:test");
const assert=require("node:assert/strict");
const {DEFAULT_JOB,usesGPU,numericInput,walltimeSeconds,jobErrors,clusterPayload,bundleFilename}=require("../static/cluster-view.js");

test("cluster defaults request independent GPU replicas with scheduler-controlled devices",()=>{
  assert.equal(DEFAULT_JOB.platform,"CUDA");
  assert.equal(DEFAULT_JOB.replicas,4);
  assert.equal(DEFAULT_JOB.device_index,"");
  assert.equal(DEFAULT_JOB.snapshot_every,1000);
  assert.deepEqual(jobErrors(DEFAULT_JOB),[]);
  assert.equal(usesGPU("OpenCL"),true);
  assert.equal(usesGPU("CPU"),false);
  assert.equal(usesGPU("numpy"),false);
});

test("blank numeric edits stay invalid instead of becoming zero or an old setting",()=>{
  assert.equal(numericInput("  "),null);
  assert.equal(numericInput(null),null);
  assert.equal(numericInput("0"),0);
  for(const key of ["steps","replicas","sample_every","snapshot_every","cpus","memory_gb"]){
    assert.ok(jobErrors({...DEFAULT_JOB,[key]:null}).length>0,key);
    assert.ok(jobErrors({...DEFAULT_JOB,[key]:""}).length>0,key);
  }
  assert.equal(jobErrors({...DEFAULT_JOB,snapshot_every:0}).length,0);
  assert.ok(jobErrors({...DEFAULT_JOB,sample_every:0}).length>0);
});

test("resource and output bounds reject fractions, negatives, infinity, and excess counts",()=>{
  assert.ok(jobErrors({...DEFAULT_JOB,replicas:1025}).length>0);
  assert.ok(jobErrors({...DEFAULT_JOB,steps:1.2}).length>0);
  assert.ok(jobErrors({...DEFAULT_JOB,steps:Infinity}).length>0);
  assert.ok(jobErrors({...DEFAULT_JOB,cpus:-1}).length>0);
  assert.ok(jobErrors({...DEFAULT_JOB,memory_gb:4097}).length>0);
  assert.equal(jobErrors({...DEFAULT_JOB,steps:10,snapshot_every:100}).length,0);
});

test("wall-time parsing accepts elapsed-hours and day formats without malformed durations",()=>{
  assert.equal(walltimeSeconds("01:00:00"),3600);
  assert.equal(walltimeSeconds("100:30:20"),361820);
  assert.equal(walltimeSeconds("2-03:04:05"),183845);
  assert.equal(walltimeSeconds("365-00:00:00"),365*86400);
  for(const invalid of ["00:00:00","1:60:00","1:00:60","1-24:00:00","366-00:00:00","1:2:3","tomorrow",""])
    assert.equal(walltimeSeconds(invalid),null,invalid);
});

test("a single logical GPU index is explicit and cannot leak into CPU execution",()=>{
  assert.equal(jobErrors({...DEFAULT_JOB,device_index:"0"}).length,0);
  assert.equal(jobErrors({...DEFAULT_JOB,device_index:""}).length,0);
  for(const device_index of ["0,1","-1","1.5","GPU-abc","65536"])
    assert.ok(jobErrors({...DEFAULT_JOB,device_index}).length>0,device_index);
  assert.ok(jobErrors({...DEFAULT_JOB,platform:"CPU",device_index:"0"}).length>0);
  assert.equal(jobErrors({...DEFAULT_JOB,platform:"CPU",device_index:""}).length,0);
});

test("export payload snapshots the draft without mutating it or adding live state",()=>{
  const config={name:"Unapplied model",chromatin:{beads:400},genes:[{id:"gene_a"}]};
  const job={...DEFAULT_JOB,steps:"20",account:"  lab  "};
  const payload=clusterPayload(config,job);
  assert.equal(payload.job.steps,20);
  assert.equal(payload.job.account,"lab");
  assert.equal(job.steps,"20");
  assert.equal(job.account,"  lab  ");
  assert.deepEqual(Object.keys(payload).sort(),["config","job"]);
  payload.config.genes[0].id="changed";
  payload.config.chromatin.beads=100;
  assert.equal(config.genes[0].id,"gene_a");
  assert.equal(config.chromatin.beads,400);
  assert.throws(()=>clusterPayload(config,{...job,steps:null}),/Steps per replica/);
  assert.throws(()=>clusterPayload(null,job),/draft model/);
});

test("download names honor a safe ZIP response and reject paths or non-ZIP names",()=>{
  assert.equal(bundleFilename('attachment; filename="job-42.zip"',"model"),"job-42.zip");
  assert.equal(bundleFilename('attachment; filename="../../other.zip"',"my model"),"my-model-cluster.zip");
  assert.equal(bundleFilename('attachment; filename="model.html"',"model"),"model-cluster.zip");
  assert.equal(bundleFilename(null,""),"chromatin-studio-cluster.zip");
});
