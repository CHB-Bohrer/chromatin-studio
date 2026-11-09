"use strict";
const test=require("node:test"),assert=require("node:assert/strict");
const {parseMSD,parseCurveInput,calibrationPayload,physicalPreview,syntheticExample,formatNumber,benchmarkChart}=require("../static/calibration-view.js");
const csv="lag,msd\n.1,.2\n.2,.3\n.5,.8\n1,1.2";
test("MSD input accepts explicit headers, comments, whitespace and scientific notation",()=>{
  assert.deepEqual(parseMSD("# measured\nlag,msd\n.1,2e-1\n.2,.3\n.5,.8\n1,1.2"),{lag:[.1,.2,.5,1],msd:[.2,.3,.8,1.2]});
  assert.deepEqual(parseMSD("0.1 0.2\n0.2\t0.3\n0.5;0.8\n1 1.2").lag,[.1,.2,.5,1]);
});
test("MSD parser rejects zero log coordinates, duplicate lags and extra columns",()=>{
  assert.throws(()=>parseMSD("0,0\n.1,.2\n.2,.3\n.5,.8"),/positive finite/);
  assert.throws(()=>parseMSD(".1,.2\n.1,.3\n.5,.8\n1,1.2"),/strictly increasing/);
  assert.throws(()=>parseMSD(".1,.2,5\n.2,.3\n.5,.8\n1,1.2"),/exactly two/);
  assert.throws(()=>parseMSD(".1,.2\n.2,.3\n.5,.8"),/at least four/);
  assert.throws(()=>parseMSD(Array.from({length:513},(_,i)=>(i+1)+",1").join("\n")),/at most 512/);
});
test("CLI JSON preserves measurement conventions and baseline provenance",()=>{
  const raw={lag:[.1,.2,.5,1],msd:[.2,.3,.8,1.2],dimension:2,observable:"relative_pair",drift_correction:"none",source:"run-17",baseline_fingerprint:"a".repeat(64),counts:[20,18,14,10],lag_unit:"reduced_time",msd_unit:"reduced_length_squared"};
  assert.deepEqual(parseCurveInput(JSON.stringify(raw)),raw);
  assert.equal(parseCurveInput(JSON.stringify({curve:raw})).observable,"relative_pair");
});
test("JSON also accepts simple coordinate rows",()=>{
  assert.deepEqual(parseCurveInput(JSON.stringify([[.1,.2],[.2,.3],[.5,.8],[1,1.2]])).lag,[.1,.2,.5,1]);
  assert.throws(()=>parseCurveInput('{"lag":[1,2],"msd":[1]}'),/matching/);
});
test("payload checks JSON conventions instead of silently relabeling a curve",()=>{
  const f=syntheticExample(),curve={...parseMSD(csv),dimension:2,observable:"locus",drift_correction:"none"};
  assert.throws(()=>calibrationPayload({...f,simulation:JSON.stringify(curve)}),/conflicts/);
  const valid=calibrationPayload({...f,dimension:"2",simulation:JSON.stringify(curve)});
  assert.equal(valid.simulation.dimension,2);
  assert.throws(()=>calibrationPayload({...f,observable:"relative_pair",drift_correction:"center_of_mass"}),/already cancels/);
});
test("physical input units and positive independently sourced length are required",()=>{
  const f=syntheticExample();
  assert.throws(()=>calibrationPayload({...f,length_nm:""}),/length/);
  assert.throws(()=>calibrationPayload({...f,length_source:""}),/independent evidence/);
  assert.throws(()=>calibrationPayload({...f,experiment_source:""}),/source/);
  assert.throws(()=>calibrationPayload({...f,fit_min_s:"20",fit_max_s:"1"}),/increasing/);
  assert.throws(()=>calibrationPayload({...f,simulation:JSON.stringify({...parseMSD(csv),lag_unit:"s"})}),/incompatible/);
});
test("payload snapshots draft state and always retains synthetic classification",()=>{
  const config={name:"draft",chromatin:{beads:12}},f=syntheticExample(),p=calibrationPayload(f,config);
  assert.equal(p.source_kind,"synthetic");p.config.chromatin.beads=100;assert.equal(config.chromatin.beads,12);
  assert.equal(p.fit_min_s,1);assert.equal(p.fit_max_s,20);
});
test("JSON source metadata can supply provenance without discarding its fingerprint",()=>{
  const f=syntheticExample(),curve={...parseMSD(csv),source:"trajectory/run-1",baseline_fingerprint:"b".repeat(64)};
  const p=calibrationPayload({...f,simulation:JSON.stringify(curve),simulation_source:""});
  assert.equal(p.simulation.source,"trajectory/run-1");assert.equal(p.simulation.baseline_fingerprint,"b".repeat(64));
});
test("physical previews apply one clock to step, diffusion and elongation",()=>{
  const c={simulation:{dt:.002},chromatin:{diffusion:.2,bp_per_bead:200},proteins:[{id:"tf",diffusion:1}],genes:[{id:"g",elongation_rate:2}]};
  const rows=physicalPreview(c,100,2);
  assert.equal(rows[0].value,.004);assert.ok(Math.abs(rows[1].value-.001)<1e-14);
  assert.ok(Math.abs(rows[2].value-.005)<1e-14);assert.equal(rows[3].value,200);
  assert.deepEqual(physicalPreview(c,100,0),[]);
});
test("benchmark pass counts and very small errors format without precision exceptions",()=>{
  assert.equal(formatNumber(14,0),"14");
  assert.match(formatNumber(1e-11),/e-11/);
  assert.equal(formatNumber(Infinity),"—");
});
test("saved benchmark curves render sampled and analytical values as finite SVG coordinates",()=>{
  const r=require("../benchmarks/reference-results.json"),html=benchmarkChart(r.curves);
  assert.ok(html.includes("Rouse"));assert.ok(html.includes('stroke-dasharray="5 4"'));
  assert.ok(html.includes('role="img"'));assert.ok(!html.includes("NaN"));assert.ok(!html.includes("Infinity"));
});
test("public reference import removes synthetic simulation and length evidence",()=>{
  const {referenceFields}=require('../static/calibration-view.js');
  const experiment={lag:[20,40,60,80],msd:[.1,.2,.3,.4],lag_unit:'s',msd_unit:'um2',dimension:3,observable:'relative_pair',drift_correction:'none',source:'public archive',analysis:{condition:'RAD21 2h',partition:'fit'},counts:[10,9,8,7]};
  const fields=referenceFields({reference:{id:'reference',title:'Measured',default_fit_min_s:20,default_fit_max_s:80},experiment},syntheticExample());
  assert.equal(fields.source_kind,'experimental');assert.equal(fields.simulation,'');assert.equal(fields.length_nm,'');assert.equal(fields.simulation_source,'');assert.equal(fields.length_source,'');
  assert.equal(fields.observable,'relative_pair');assert.deepEqual(JSON.parse(fields.experiment),experiment);
  assert.throws(()=>calibrationPayload(fields,{}),/Simulation/);
});
test("public reference retains real simulation provenance and rejects missing unit declarations",()=>{
  const {referenceFields}=require('../static/calibration-view.js');
  const experiment={lag:[20,40,60,80],msd:[.1,.2,.3,.4],lag_unit:'s',msd_unit:'um2',dimension:3,observable:'relative_pair',drift_correction:'none'};
  const data={reference:{id:'heldout',title:'Held out'},experiment};
  const f=referenceFields(data,{simulation:'real curve',simulation_source:'run-1',length_nm:'100',source_kind:'experimental'});
  assert.equal(f.simulation_source,'run-1');assert.equal(f.length_nm,'100');
  assert.throws(()=>referenceFields({...data,experiment:{...experiment,lag_unit:'reduced_time'}}),/seconds/);
});
test('measured reference plot displays actual data without claiming a fitted clock',()=>{
  const {referenceChart}=require('../static/calibration-view.js');
  const c=require('../references/hansen/curves/fbn2-rad21-2h-fit-corrected.json');
  const html=referenceChart(c);
  assert.ok(html.includes('Experimental curve only'));
  assert.ok(html.includes('no simulated curve or time mapping'));
  assert.ok(html.includes('seconds (log scale)'));
  assert.ok(!html.includes('NaN'));assert.ok(!html.includes('Infinity'));
});
