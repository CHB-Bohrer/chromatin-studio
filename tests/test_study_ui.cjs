"use strict";
const test=require("node:test"),assert=require("node:assert/strict"),fs=require("node:fs"),vm=require("node:vm");
const {finite,metricValue,compatibleConditions,conditionMetrics,checkedOptions,effectExtent}=require("../static/study-view.js");
const catalog={engines:[{id:"reference",name:"Four-promoter reference"},{id:"polymer",name:"Polymer pilot"}],conditions:[{id:"baseline",name:"Baseline",engines:["reference","polymer"]},{id:"uncoupled",name:"Independent promoters",engines:["reference"]},{id:"no_cohesin",name:"No cohesin",engines:["polymer"]}],defaults:{reference:{engine:"reference",conditions:["baseline","uncoupled"],cells:64,duration:2,seed:42},polymer:{engine:"polymer",conditions:["baseline","no_cohesin"],replicas:2}},limits:{reference:{max_work:200000},polymer:{max_particle_steps:20000000}}};
function browserHarness({session,fetch=async()=>{throw Error("Unexpected network request");},config={name:"Original draft",genes:[{id:"a"},{id:"b"}]}}={}){
  const root={innerHTML:"",listeners:{},addEventListener(name,fn){this.listeners[name]=fn;},removeEventListener(name){delete this.listeners[name];},querySelector(){return null;}},window={},navigation=[];
  vm.runInNewContext(fs.readFileSync(require.resolve("../static/study-view.js"),"utf8"),{window,document:{},fetch,setTimeout,Blob,URL});
  const state=session||{catalog,options:{...catalog.defaults.reference}};
  const view=window.ChromatinStudy.render(root,{getConfig:()=>config,session:state,onNavigate:target=>navigation.push(target)});
  return {root,view,session:state,config,navigation};
}
test("missing and invalid measurements remain undefined while zero remains a real value",()=>{
  for(const value of [null,undefined,NaN,Infinity,"0"]){assert.equal(finite(value),false);assert.equal(metricValue(value),"Undefined");}
  assert.equal(metricValue(0),"0");assert.equal(metricValue(-.125),"-0.125");
  const m=conditionMetrics({metrics:{observed:{phi_cis:0,phi_trans:null,delta_phi:null,reason:"Constant allele"}}},"reference");
  assert.equal(m.cis,0);assert.equal(m.delta,null);assert.equal(m.trans,null);assert.equal(m.reason,"Constant allele");
});
test("reference controls stay optional while polymer baseline is included without mutating options",()=>{
  assert.deepEqual(compatibleConditions(catalog,"polymer").map(c=>c.id),["baseline","no_cohesin"]);
  const options={engine:"reference",conditions:["uncoupled","uncoupled"]};
  assert.deepEqual(checkedOptions(catalog,options).conditions,["uncoupled"]);
  assert.deepEqual(checkedOptions(catalog,{engine:"polymer",conditions:["no_cohesin"]}).conditions,["baseline","no_cohesin"]);
  assert.throws(()=>checkedOptions(catalog,{engine:"reference",conditions:[]}),/at least one/);
  assert.deepEqual(options.conditions,["uncoupled","uncoupled"]);
  assert.throws(()=>checkedOptions(catalog,{engine:"reference",conditions:["no_cohesin"]}),/not available/);
});
test("polymer measurements never fabricate a matched trans contrast",()=>{
  const m=conditionMetrics({metrics:{nascent_phi:.4,promoter_phi:.1,count_correlation:.3,phi_trans:.8,delta_phi:-.4,bootstrap_ci95:[.1,.6]}},"polymer");
  assert.equal(m.nascent,.4);assert.equal(m.promoter,.1);assert.equal(m.cis,null);assert.equal(m.trans,null);assert.equal(m.delta,null);assert.deepEqual(m.interval,[.1,.6]);
});
test("contrast plot includes a delta above one and its uncertainty instead of clipping",()=>{
  const conditions=[{metrics:{observed:{phi_cis:1,phi_trans:-1,delta_phi:2,ci95:[1.2,2]}}}];
  assert.equal(effectExtent(conditions),2);assert.equal(effectExtent([]),.25);
  assert.equal(effectExtent([{metrics:{observed:{delta_phi:-.28,ci95:[-.7,null]}}}]),.75);
});
test("rendered results expose uncertainty, condition-level warnings and undefined correlations",()=>{
  const session={catalog,options:{...catalog.defaults.reference},layer:"observed",document:{kind:"run",configName:"Draft <script>",data:{engine:"reference",options:{cells:64},conditions:[{id:"baseline",name:"Baseline",seed:42,warnings:["RNA capacity blocked initiation"],metrics:{observed:{phi_cis:null,phi_trans:0,delta_phi:null,ci95:[null,null],reason:"Constant allele"},mean_nascent_probability:[0,0,1,1],mean_nascent_counts:[0,0,2,2]}}]}}};
  const {root}=browserHarness({session});
  assert.match(root.innerHTML,/RNA capacity blocked initiation/);assert.match(root.innerHTML,/Constant allele/);assert.match(root.innerHTML,/Bootstrap support is not a p-value/);
  assert.match(root.innerHTML,/Draft &lt;script&gt;/);assert.doesNotMatch(root.innerHTML,/(?:NaN|Infinity)%/);
  assert.doesNotMatch(root.innerHTML,/<i class="study-effect-dot"/);
});
test("back navigation returns to the real system editor",()=>{
  const h=browserHarness();h.root.listeners.click({target:{closest:()=>({dataset:{studyAction:"back"}})}});
  assert.deepEqual(h.navigation,["system"]);
});
test("independent study snapshots the draft, uses only study API and survives leaving the workspace",async()=>{
  const requests=[];let finish;const deferred=new Promise(resolve=>finish=resolve);
  const h=browserHarness({fetch:async(path,options)=>{requests.push({path,body:JSON.parse(options.body)});await deferred;return {ok:true,json:async()=>({engine:"reference",options:{cells:64},conditions:[]})};}});
  const pending=h.view.execute("run");assert.equal(h.session.busy,"run");h.config.name="Later draft";h.view.destroy();finish();await pending;
  assert.equal(requests.length,1);assert.equal(requests[0].path,"/api/studies/run");assert.equal(requests[0].body.config.name,"Original draft");
  assert.equal(h.config.name,"Later draft");assert.equal(h.session.document.configName,"Original draft");assert.equal(h.session.busy,null);assert.equal(h.root.innerHTML,"");
  const resumed=browserHarness({session:h.session,config:h.config});assert.match(resumed.root.innerHTML,/Original draft/);
});
test("failed pilot keeps the prior results available and reports the backend reason",async()=>{
  const previous={kind:"plan",configName:"Previous",data:{engine:"reference",options:{},conditions:[]}};
  const session={catalog,options:{...catalog.defaults.reference},document:previous};
  const h=browserHarness({session,fetch:async()=>({ok:false,json:async()=>({error:"Pilot budget exceeded",details:["Reduce cells"]})})});
  await h.view.execute("run");assert.equal(session.document,previous);assert.equal(session.busy,null);assert.match(h.root.innerHTML,/Pilot budget exceeded/);assert.match(h.root.innerHTML,/Reduce cells/);assert.match(h.root.innerHTML,/Previous/);
});


test("reference baseline can be deselected and the request includes only the chosen controls",async()=>{
  let payload;
  const h=browserHarness({fetch:async(path,options)=>{payload=JSON.parse(options.body);return {ok:true,json:async()=>({engine:"reference",conditions:[],options:payload.options})};}});
  const baselineInput=h.root.innerHTML.match(/<input[^>]+data-study-condition="baseline"[^>]*>/)?.[0];
  assert.ok(baselineInput);assert.doesNotMatch(baselineInput,/disabled/);
  h.view.changeInput({target:{dataset:{studyCondition:"baseline"},checked:false}});
  await h.view.execute("run");assert.deepEqual(payload.options.conditions,["uncoupled"]);
  h.view.changeInput({target:{dataset:{studyCondition:"uncoupled"},checked:false}});
  await h.view.execute("run");assert.match(h.session.error,/at least one/);assert.deepEqual(payload.options.conditions,["uncoupled"]);
});


test("typed numeric options are submitted immediately without a blur or focused-field rerender",async()=>{
  let payload;
  const h=browserHarness({fetch:async(path,options)=>{payload=JSON.parse(options.body);return {ok:true,json:async()=>({engine:"reference",conditions:[],options:payload.options})};}});
  const before=h.root.innerHTML;
  h.root.listeners.input({target:{type:"number",tagName:"INPUT",dataset:{studyOption:"duration"},value:"24"}});
  assert.equal(h.session.options.duration,24);assert.equal(h.root.innerHTML,before);
  await h.view.execute("plan");assert.equal(payload.options.duration,24);
  h.root.listeners.input({target:{type:"number",tagName:"INPUT",dataset:{studyOption:"duration"},value:""}});
  assert.equal(h.session.options.duration,null);
  h.view.destroy();assert.equal(h.root.listeners.input,undefined);
});
