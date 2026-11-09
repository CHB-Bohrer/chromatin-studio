"use strict";
const test=require("node:test"),assert=require("node:assert/strict");
const {wallAt,bodyFactor,transcriptionObservation}=require("../static/model-explain.js");
test("soft wall is force-free inside and its outside force is the energy derivative",()=>{
  assert.deepEqual(wallAt(10,30,8),{penetration:0,energy:0,inwardForce:0});
  assert.deepEqual(wallAt(10,30,10),{penetration:0,energy:0,inwardForce:0});
  assert.deepEqual(wallAt(10,30,12),{penetration:2,energy:60,inwardForce:60});
  const eps=1e-5,derivative=(wallAt(10,30,12+eps).energy-wallAt(10,30,12-eps).energy)/(2*eps);
  assert.ok(Math.abs(derivative-60)<1e-7);
  assert.deepEqual(wallAt(10,0,12),{penetration:2,energy:0,inwardForce:0});
  assert.equal(wallAt(0,30,2),null);assert.equal(wallAt(10,-1,2),null);assert.equal(wallAt(10,null,2),null);
});
test("body response has the requested center value, neutral far field, and disabled-field behavior",()=>{
  assert.equal(bodyFactor(0,2,.25).factor,.25);assert.equal(bodyFactor(0,2,4).factor,4);
  const atRadius=bodyFactor(2,2,4);assert.ok(Math.abs(atRadius.weight-Math.exp(-.5))<1e-15);
  assert.ok(Math.abs(Math.log(atRadius.factor)-Math.exp(-.5)*Math.log(4))<1e-15);
  assert.equal(bodyFactor(100,2,4).factor,1);assert.equal(bodyFactor(2,2,1).factor,1);
  assert.deepEqual(bodyFactor(0,2,4,false),{weight:0,factor:1});
  assert.equal(bodyFactor(0,0,2),null);assert.equal(bodyFactor(0,2,0),null);assert.equal(bodyFactor(0,2,NaN),null);
});
test("promoter state and nascent positivity are distinct observations, including missing values",()=>{
  assert.deepEqual(transcriptionObservation({on:false,nascent:1}),{promoterOn:false,nascentPositive:true});
  assert.deepEqual(transcriptionObservation({on:true,nascent:0}),{promoterOn:true,nascentPositive:false});
  assert.deepEqual(transcriptionObservation({on:false,nascent:null}),{promoterOn:false,nascentPositive:null});
  assert.deepEqual(transcriptionObservation(null),{promoterOn:null,nascentPositive:null});
});
