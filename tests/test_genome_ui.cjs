"use strict";

const assert = require("node:assert/strict");
const test = require("node:test");
const genome = require("../static/genome-view.js");

test("20 Mb selection respects bin boundaries and clamps the final axis endpoint", () => {
  const bead = bp => genome.bpToBead(bp, 10_000, 2000);
  assert.equal(bead(0), 0);
  assert.equal(bead(9_999), 0);
  assert.equal(bead(10_000), 1);
  assert.equal(bead(19_999_999), 1999);
  assert.equal(bead(20_000_000), 1999);
  assert.equal(bead(-1), 0);
  assert.equal(bead(20_000_001), 1999);
});

test("inclusive gene bead endpoints describe the complete half-open genomic interval", () => {
  assert.deepEqual(genome.beadInterval(0, 10_000), {start: 0, end: 10_000});
  assert.deepEqual(genome.beadInterval(1999, 10_000), {start: 19_990_000, end: 20_000_000});
  // The preset's gene A includes both beads 330 and 341: twelve bins, 120 kb.
  const first = genome.beadInterval(330, 10_000);
  const last = genome.beadInterval(341, 10_000);
  assert.equal(first.start, 3_300_000);
  assert.equal(last.end, 3_420_000);
  assert.equal(last.end - first.start, 120_000);
});

test("RLE boundary lookup never assigns the previous run to the next bead", () => {
  const runs = [[0, 2, 0], [2, 3, 4], [3, 7, 0.25], [7, 8, 2]];
  assert.deepEqual(Array.from({length: 8}, (_, i) => genome.rleValue(runs, i)),
    [0, 0, 4, 0.25, 0.25, 0.25, 0.25, 2]);
  assert.equal(genome.rleValue(runs, -1), null);
  assert.equal(genome.rleValue(runs, 8), null);
  assert.equal(genome.rleValue([], 0), null);
  assert.equal(genome.rleValue([[0, 1, 3], [2, 3, 4]], 1), null);
});

test("visible runs clip at both window edges without changing authoritative data", () => {
  const runs = [[0, 2, 0], [2, 3, 4], [3, 7, 0.25], [7, 8, 2]];
  const original = JSON.stringify(runs);
  assert.deepEqual(genome.visibleRuns(runs, 1, 5), [[1, 2, 0], [2, 3, 4], [3, 5, 0.25]]);
  assert.deepEqual(genome.visibleRuns(runs, 2, 3), [[2, 3, 4]]);
  assert.deepEqual(genome.visibleRuns(runs, 7, 8), [[7, 8, 2]]);
  assert.deepEqual(genome.visibleRuns(runs, 8, 9), []);
  assert.deepEqual(genome.visibleRuns(runs, 3, 3), []);
  assert.deepEqual(genome.visibleRuns([], 0, 8), []);
  assert.equal(JSON.stringify(runs), original);
});

test("panning preserves window width while staying within genomic bounds", () => {
  assert.deepEqual(genome.clampWindow(100, 500, 2000, 100), {start: 100, end: 500});
  assert.deepEqual(genome.clampWindow(-100, 300, 2000, 100), {start: 0, end: 400});
  assert.deepEqual(genome.clampWindow(1800, 2200, 2000, 100), {start: 1600, end: 2000});
  assert.deepEqual(genome.clampWindow(-1000, 3000, 2000, 100), {start: 0, end: 2000});
});

test("zoomed windows remain selectable and cannot be narrower than the requested minimum", () => {
  for (const [start, end] of [[0, 0], [1999, 2000], [750, 751], [1900, 1850]]) {
    const window = genome.clampWindow(start, end, 2000, 100);
    assert.ok(window.start >= 0);
    assert.ok(window.end <= 2000);
    assert.ok(window.end - window.start >= 100);
  }
  assert.deepEqual(genome.clampWindow(0, 1, 50, 100), {start: 0, end: 50});
});

test("local occupancy and persistence summaries weight clipped genomic lengths", () => {
  assert.deepEqual(genome.summarizeRuns([[0, 2, 0], [2, 6, 1]], 1, 5),
    {count: 4, mean: 0.75, min: 0, max: 1});
  assert.deepEqual(genome.summarizeRuns([[0, 2, 1], [2, 6, 3]], 1, 5),
    {count: 4, mean: 2.5, min: 1, max: 3});
  assert.deepEqual(genome.summarizeRuns([[0, 2, 0]], 0, 1),
    {count: 1, mean: 0, min: 0, max: 0});
  assert.equal(genome.summarizeRuns([], 0, 2), null);
  assert.equal(genome.summarizeRuns([[0, 2, 1]], 2, 4), null);
});

test("nucleosome glyphs require both nucleosome-scale resolution and enough screen space", () => {
  assert.equal(genome.nucleosomeDisplay(200, 24), "individual");
  assert.equal(genome.nucleosomeDisplay(200, 8), "aggregate");
  assert.equal(genome.nucleosomeDisplay(10_000, 80), "coarse");
  assert.equal(genome.nucleosomeDisplay(10_000, 2), "coarse");
});

test("direct molecule glyphs deduplicate IDR domains while retaining separate binding sites", () => {
  const sites = [{bead: 2, molecule: 7, species: "tf", particle: 40},
    {bead: 2, molecule: 7, species: "tf", particle: 41},
    {bead: 3, molecule: 7, species: "tf", particle: 42},
    {bead: 3, molecule: 8, species: "tf", particle: 44},
    {bead: 3, molecule: 9, species: "remodeler", particle: 46}];
  assert.deepEqual(genome.visibleBindingSites(sites, 2, 4, "tf"), [
    {bead: 2, molecule: 7, species: "tf", particles: [40, 41]},
    {bead: 3, molecule: 7, species: "tf", particles: [42]},
    {bead: 3, molecule: 8, species: "tf", particles: [44]}]);
  assert.equal(genome.visibleBindingSites(sites, 3, 4).length, 3);
  assert.deepEqual(genome.visibleBindingSites(sites, 4, 5), []);
});

test("recorded events retain engine IDs and observation order, never infer missed transitions", () => {
  const stream = {run_id: "a", first_retained_id: 1, last_id: 2, events: [
    {id: 1, time: 1.0, type: "chromatin_bind"},
    {id: 2, time: 0.9, type: "nucleosome_eviction"}]};
  let state = genome.mergeRecordedEvents(null, stream);
  assert.deepEqual(state.events.map(e => e.id), [1, 2]);
  assert.deepEqual(state.events.map(e => e.time), [1, .9]);
  assert.equal(state.missing, 0);
  state = genome.mergeRecordedEvents(state, stream);
  assert.equal(state.events.length, 2);
  state = genome.mergeRecordedEvents(state, {run_id:"a", first_retained_id:6,last_id:7,
    events:[{id:6,type:"chromatin_unbind"},{id:7,type:"chromatin_bind"}]});
  assert.deepEqual(state.events.map(e => e.id), [1, 2, 6, 7]);
  assert.equal(state.missing, 3);
  state = genome.mergeRecordedEvents(state, {run_id:"a",first_retained_id:6,last_id:7,
    events:[{id:6},{id:7}]});
  assert.equal(state.missing, 3, "repeated snapshots must not count lost events again");
});

test("event buffers reset between applied simulations even when time is unchanged", () => {
  const old = genome.mergeRecordedEvents(null, {run_id:"old",last_id:5,
    first_retained_id:4,events:[{id:4},{id:5}]});
  assert.equal(old.missing, 3);
  const fresh = genome.mergeRecordedEvents(old,{run_id:"new",last_id:1,
    first_retained_id:1,events:[{id:1,time:0}]});
  assert.deepEqual(fresh.events,[{id:1,time:0}]);
  assert.equal(fresh.missing,0);
  assert.equal(genome.mergeRecordedEvents(fresh,undefined).available,false);
});

test("display history is bounded without labelling intentionally pruned events missing", () => {
  const state=genome.mergeRecordedEvents(null,{run_id:"a",first_retained_id:1,last_id:300,
    events:Array.from({length:300},(_,i)=>({id:i+1}))},240);
  assert.equal(state.events.length,240);
  assert.equal(state.events[0].id,61);
  assert.equal(state.missing,0);
});

test("event filters include a nucleosome sliding into the view and preserve species-free turnover", () => {
  const events=[{id:1,bead:2,to_bead:3,species:"remodeler"},
    {id:2,bead:3,species:"tf"},{id:3,bead:3},{id:4,bead:5,species:"tf"}];
  assert.deepEqual(genome.filterRecordedEvents(events,3,4).map(e=>e.id),[1,2,3]);
  assert.deepEqual(genome.filterRecordedEvents(events,3,4,new Set(["tf"])).map(e=>e.id),[2,3]);
});

test("selected site uses exact current kinetic tracks including zero hazard", () => {
  const protein={id:"tf",on_rate:2,off_rate:.5};
  const tracks={on_rate_runs:{tf:[[0,1,0],[1,2,8]]},off_rate_runs:{tf:[[0,1,2],[1,2,.2]]},affinity_runs:{tf:[[0,2,4]]}};
  assert.equal(genome.currentRate(tracks,protein,0,"on"),0);
  assert.equal(genome.currentRate(tracks,protein,1,"on"),8);
  assert.equal(genome.currentRate(tracks,protein,1,"off"),.2);
  assert.equal(genome.currentRate({affinity_runs:tracks.affinity_runs},protein,0,"on"),8);
  assert.equal(genome.currentRate({affinity_runs:tracks.affinity_runs},protein,0,"off"),null,
    "an old snapshot cannot establish occupancy-dependent release hazards");
});

function liveGenomeHarness(){
  const vm=require("node:vm"),fs=require("node:fs"),nodes=new Map(),drawn={arcs:[],text:[]},callbacks={run:0,step:0};
  const context=new Proxy({measureText:text=>({width:String(text).length*5}),arc:(...args)=>drawn.arcs.push(args),fillText:text=>drawn.text.push(String(text))},{get:(target,key)=>target[key]||(()=>{})});
  function node(key){if(!nodes.has(key))nodes.set(key,{innerHTML:"",textContent:"",value:"",checked:true,disabled:false,hidden:false,dataset:{},style:{},listeners:{},attributes:{},classList:{toggle(){}},addEventListener(name,fn){this.listeners[name]=fn;},setAttribute(name,value){this.attributes[name]=value;},querySelectorAll:()=>[],setCustomValidity(){},reportValidity(){}});return nodes.get(key);}
  const canvas=node('[data-role="canvas"]');Object.assign(canvas,{parentElement:{clientWidth:900},getContext:()=>context});
  const root={innerHTML:"",hidden:false,listeners:{},classList:{toggle(){}},addEventListener(name,fn){this.listeners[name]=fn;},querySelector:key=>node(key)};
  const window={devicePixelRatio:1};
  vm.runInNewContext(fs.readFileSync(require.resolve("../static/genome-view.js"),"utf8"),{window,document:{activeElement:null},ResizeObserver:class{observe(){}},requestAnimationFrame:()=>{}});
  const view=new window.ChromatinGenomeView(root,{onRun:()=>callbacks.run++,onStep:()=>callbacks.step++});
  const config={chromatin:{beads:8,bp_per_bead:200},genes:[],regulatory_regions:[],marks:[],proteins:[{id:"tf",name:"TF <safe>",color:"#abcdef",on_rate:2,off_rate:.2}]};
  const snapshot={time:1,step:500,genes:[],nucleosomes:{enabled:true},genome_tracks:{nucleosome_occupancy_runs:[[0,3,1],[3,5,0],[5,8,1]],persistence_runs:[[0,8,2]],affinity_runs:{tf:[[0,3,.5],[3,5,2],[5,8,.5]]},occupancy_runs:{tf:[[0,3,0],[3,4,1],[4,8,0]]},binding_sites:[{bead:3,molecule:9,particle:20,species:"tf",off_rate:.1}],on_rate_runs:{tf:[[0,3,1],[3,5,4],[5,8,1]]},off_rate_runs:{tf:[[0,3,.5],[3,5,.1],[5,8,.5]]}},interaction_state:{species:[{id:"tf",direct_bound:1,free:9}]},genomic_events:{run_id:"test",first_retained_id:1,last_id:3,events:[{id:1,type:"chromatin_bind",species:"tf",molecule:9,bead:3,time:.7},{id:2,type:"chromatin_unbind",species:"tf",molecule:8,bead:4,time:.8},{id:3,type:"nucleosome_reset",bead:5,time:1,interval_start:.98,interval_end:1,timing:"interval_endpoint"}]}};
  return {view,root,node,drawn,callbacks,config,snapshot};
}

test("live controls remain paused until the user acts and follow applied run availability",()=>{
  const h=liveGenomeHarness();h.view.update(h.snapshot,h.config,{type:"chromatin",index:3},{ready:true});
  assert.equal(h.callbacks.run,0);assert.equal(h.callbacks.step,0);
  assert.equal(h.node('[data-action="live-run"]').disabled,false);
  h.view.action("live-run");assert.equal(h.callbacks.run,1);
  h.view.action("live-step");assert.equal(h.callbacks.step,1);
  h.view.setRunState({running:true,busy:true,ready:true,dirty:false});
  assert.equal(h.node('[data-action="live-run"]').textContent,"Pause");
  assert.equal(h.node('[data-action="live-run"]').disabled,false);
  assert.equal(h.node('[data-action="live-step"]').disabled,true);
  h.view.setRunState({running:false,busy:false,dirty:true});
  assert.equal(h.node('[data-action="live-run"]').disabled,true);
});

test("accessible event feed and selected-site rates use live records with interval and release explanations",()=>{
  const h=liveGenomeHarness();h.view.update(h.snapshot,h.config,{type:"chromatin",index:3},{ready:true});
  assert.match(h.node('[data-role="event-status"]').textContent,/Latest 3 events cached \(limit 240\); 3 shown/);
  assert.match(h.node('[data-role="event-status"]').textContent,/Earlier cached events are discarded/);
  const feed=h.node('[data-role="event-list"]').innerHTML;
  assert.match(feed,/TF &lt;safe&gt;/);assert.match(feed,/Released into 3D; shown at its last binding site/);
  assert.match(feed,/Interval 0.98–1 t · endpoint observation/);
  assert.ok(feed.indexOf('data-genomic-event="3"')<feed.indexOf('data-genomic-event="1"'));
  const table=h.node('[data-role="selection"]').innerHTML;
  assert.match(table,/Vacant/);assert.match(table,/<td>4<\/td><td>0.1<\/td>/);
  assert.match(table,/conditional hazard/);assert.doesNotMatch(table,/NaN|Infinity/);
});

test("live drawing retains molecule ID labels and explicitly downgrades coarse bins",()=>{
  const h=liveGenomeHarness();h.view.update(h.snapshot,h.config,null,{ready:true});h.view.draw();
  assert.ok(h.drawn.text.includes("#9"));
  assert.match(h.node('[data-role="resolution"]').innerHTML,/One circle = one modeled nucleosome site/);
  h.config.chromatin.bp_per_bead=10_000;h.view.update(h.snapshot,h.config,null);h.view.draw();
  assert.match(h.node('[data-role="resolution"]').innerHTML,/individual nucleosomes are unresolved/);
  h.view.action("all-tracks");assert.equal(h.view.live,false);
  h.view.draw();assert.ok(h.drawn.text.includes("Persistence"));
});
