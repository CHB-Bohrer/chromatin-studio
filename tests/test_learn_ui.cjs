"use strict";
const test=require("node:test"),assert=require("node:assert/strict");
const {clamp,diffusionMoments,eventProbability,idealBinding,chromatinSite,regulatoryRates,rnaProgress}=require("../static/learn-view.js");
test("Brownian illustration separates thermal MSD from force-induced drift",()=>{
  const free=diffusionMoments(.5,2),forced=diffusionMoments(.5,2,-3);
  assert.equal(free.msd,6);assert.equal(free.axisVariance,2);assert.equal(free.mean,0);
  assert.equal(forced.msd,free.msd);assert.equal(forced.mean,-3);
  assert.deepEqual(diffusionMoments(0,10,2),{msd:0,axisVariance:0,rms:0,mean:0});
});
test("finite-interval reaction probability uses exponential survival including tiny rates",()=>{
  assert.equal(eventProbability(0,2),0);
  assert.ok(Math.abs(eventProbability(2,.5)-(1-Math.exp(-1)))<1e-15);
  assert.ok(Math.abs(eventProbability(1e-20,1)/1e-20-1)<1e-14);
  assert.equal(eventProbability(1e9,1),1);
});
test("ideal unlike-pair estimate converts nm and micromolar to molar kinetic rates",()=>{
  const r=idealBinding(1,5,2),expected=4*Math.PI*.005*6.02214076e8;
  assert.ok(Math.abs(r.kon/expected-1)<1e-14);
  assert.equal(r.koff,r.kon*2e-6);assert.equal(r.residence,1/r.koff);
  assert.equal(idealBinding(1,5,2,.25).kon,r.kon*.25);
  assert.equal(idealBinding(0,5,2).residence,null);
});
test("occupancy capture factors and persistence are independent; transcription composes",()=>{
  const p={basePersistence:4,occupiedPersistence:2,vacantPersistence:.25,baseAffinity:3,occupiedAffinity:.1,vacantAffinity:1,onRate:2,transcriptionMultiplier:.5};
  const occupied=chromatinSite({...p,occupied:true}),vacant=chromatinSite({...p,occupied:false,transcribing:true});
  assert.equal(occupied.persistence,8);assert.ok(Math.abs(occupied.hazard-.6)<1e-14);
  assert.equal(vacant.persistence,.5);assert.equal(vacant.hazard,6);
});
test("lost contact removes regional gains; mark gate affects activation",()=>{
  const p={tf:2,coactivators:1,coactivatorGain:3,activationGain:.4,initiationGain:.5,onRate:.2,initiationRate:2};
  const r=regulatoryRates(p);assert.equal(r.signal,5);assert.ok(Math.abs(r.onRate-.6)<1e-14);assert.equal(r.initiationRate,7);
  const lost=regulatoryRates({...p,contact:false});assert.equal(lost.onRate,.2);assert.equal(lost.initiationRate,2);
  const gate=regulatoryRates({...p,gate:false});assert.equal(gate.onRate,0);assert.equal(gate.initiationRate,7);
});
test("RNA progress caps at full gene length and zero speed remains nascent",()=>{
  assert.deepEqual(rnaProgress(200,3,1000),{length:600,fraction:.6,released:false});
  assert.deepEqual(rnaProgress(200,20,1000),{length:1000,fraction:1,released:true});
  assert.deepEqual(rnaProgress(0,20,1000),{length:0,fraction:0,released:false});
});
test("bounded illustrative controls reject nonfinite values",()=>{
  assert.equal(clamp(Infinity,0,2,.5),.5);assert.equal(clamp(null,0,2,.5),.5);
  assert.equal(clamp(-10,0,2),0);assert.equal(clamp(20,0,2),2);
  assert.equal(eventProbability(-1,3),0);assert.equal(diffusionMoments(NaN,3).msd,0);
});
