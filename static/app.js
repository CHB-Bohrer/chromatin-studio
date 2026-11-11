/* Chromatin Studio: local, dependency-free research workbench. */
"use strict";
(() => {
  const $ = (s) => document.querySelector(s);
  const $$ = (s) => Array.from(document.querySelectorAll(s));
  const clone = (x) => JSON.parse(JSON.stringify(x));
  const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const pretty = (n, digits = 2) => Number.isFinite(Number(n)) ? Number(n).toLocaleString(undefined, {maximumFractionDigits: digits}) : "—";
  const genomicSize = bp => Math.abs(bp)>=1e6?`${pretty(bp/1e6,3)} Mb`:Math.abs(bp)>=1000?`${pretty(bp/1000,2)} kb`:`${pretty(bp,0)} bp`;
  function genomeScale(config){const beads=Number(config?.chromatin?.beads)||0,bp=Number(config?.chromatin?.bp_per_bead)||0;return {beads,bp,total:beads*bp};}
  function geneGenomicSummary(g,bp){return `${genomicSize(g.start*bp)}–${genomicSize((g.end+1)*bp)} · ${genomicSize((g.end-g.start+1)*bp)} gene · ${pretty(g.elongation_rate*bp/1000,2)} kb / reduced time elongation`;} 
  function proteinRoles(config,id){
    const regions=(config.regulatory_regions||[]).filter(r=>r.enabled!==false),roles=[];
    if(regions.some(r=>r.tf_species===id)||(config.genes||[]).some(g=>g.activator_species===id&&!regions.some(r=>r.target_gene===g.id)))roles.push("Transcription factor");
    if(regions.some(r=>r.coactivator_species===id))roles.push("Recruited coactivator");
    const protein=config.proteins?.find(p=>p.id===id);
    if(protein&&["nucleosome_eviction_rate","nucleosome_reassembly_rate","nucleosome_sliding_rate"].some(key=>protein[key]>0))roles.push("Chromatin remodeler");
    return roles;
  }
  function regulatoryModelSummary(config){
    const factors=(config.proteins||[]).filter(p=>proteinRoles(config,p.id).includes("Transcription factor"));
    const regions=(config.regulatory_regions||[]).filter(r=>r.enabled!==false);
    const tf=factors.length?factors.map(p=>`${esc(p.name)} (${pretty(p.count,0)} molecules)`).join(", "):"None configured";
    const rna=config.transcription?.enabled?"Individual polymerases and diffusing RNA":"Aggregate RNA counts; individual RNA tracking is off";
    return `<div class="info-note"><strong>Explicit TF particles:</strong> ${tf}.<br><strong>Regulatory elements:</strong> ${regions.length?`${regions.length} enabled enhancers / promoters`:"No enabled elements; gene-wide TF binding can supply the legacy gain"}.<br><strong>Transcription:</strong> stochastic promoter ON/OFF and initiation.<br><strong>RNA model:</strong> ${rna}.</div>`;
  }
  const S = {config:null, applied:null, snapshot:null, history:[], events:[], eventKeys:new Set(), capabilities:null, presets:[], tab:"system", dirty:false, running:false, busy:false, ready:false, backend:"numpy", selected:null, camera:{yaw:-0.42,pitch:0.22,zoom:1}, hitTargets:[], view:{bonds:true,proteins:true,labels:true,rna:true,occupancy:true,nucleus:true,fitNucleus:false}, presetId:"", renderPending:false};
  const colors = ["#9b79bd","#58a59a","#c7a068","#7294bd","#bc818e","#89a76f","#bf956d","#7c8abd"];
  let toastTimer, runGeneration = 0, genomeView=null, interactionView=null, regulationView=null, clusterView=null, learnView=null, calibrationView=null,studyView=null;

  async function api(path, body) {
    const options = body === undefined ? {} : {method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(body)};
    const response = await fetch(path, options);
    let result;
    try { result = await response.json(); } catch (_) { throw new Error(`The local server returned an unreadable response (${response.status}).`); }
    if (!response.ok) {
      const details = Array.isArray(result.details) ? result.details.join("\n") : result.details;
      throw new Error([result.error || `Request failed (${response.status})`, details].filter(Boolean).join("\n"));
    }
    return result;
  }
  function showError(message) { $("#error-message").textContent=message; $("#error-banner").hidden=false; }
  function toast(message) { $("#toast").textContent=message; $("#toast").hidden=false; clearTimeout(toastTimer); toastTimer=setTimeout(()=>$("#toast").hidden=true,3600); }
  function safe(action) { return async (...args)=>{try {await action(...args);} catch(e) {pause();showError(e.message || String(e));}}; }
  function normalizeConfig(c) {
    c.marks ??= []; c.proteins ??= []; c.interactions ??= []; c.bodies ??= []; c.genes ??= []; c.regulatory_regions ??= [];
    c.chromatin ??= {}; c.chromatin.initial_marks ??= []; c.chromatin.barriers ??= []; c.chromatin.state_interactions ??= []; c.chromatin.representation ??= "coarse";
    c.chromatin.nucleosomes={enabled:false,initial_occupancy:1,initial_regions:[],eviction_rate:0,reassembly_rate:0,occupied_persistence_multiplier:1,vacant_persistence_multiplier:.5,...c.chromatin.nucleosomes};
    c.transcription={enabled:false,max_transcripts:2000,rna_diffusion:.3,mature_rna_diffusion:.15,...c.transcription};
    c.simulation ??= {}; c.simulation.confinement_stiffness??=30; c.cohesin ??= {enabled:false,count:0,step_rate:0.2,off_rate:0.01,spring_k:1};
    c.proteins.forEach(p=>{p.chromatin_affinity ??= {unmodified:0};p.nucleosome_affinity={occupied:1,vacant:1,...p.nucleosome_affinity};p.nucleosome_off_rate_multiplier={occupied:1,vacant:1,...p.nucleosome_off_rate_multiplier};p.nucleosome_eviction_rate??=0;p.nucleosome_reassembly_rate??=0;p.nucleosome_sliding_rate??=0;});
    c.genes.forEach(g=>{g.nucleosome_eviction_rate??=0;g.pause_rate??=0;g.resume_rate??=1;g.occupied_elongation_multiplier??=1;});
    c.bodies.forEach(b=>{b.protein_affinity ??= {}; b.center ??= [0,0,0];});
    return c;
  }
  function setPath(path, value) { let parent=S.config; for(const key of path.slice(0,-1)){parent[key]??={};parent=parent[key];} if(value===null&&["chromatin_affinity","protein_affinity"].includes(path.at(-2)))delete parent[path.at(-1)];else parent[path.at(-1)]=value; }
  function getPath(path) { return path.reduce((a,b)=>a?.[b],S.config); }
  function pathAttr(path) { return esc(JSON.stringify(path)); }
  function updateBodyResponsePreviews() { document.querySelectorAll("[data-body-response-id]").forEach(row=>{const body=S.config.bodies?.find(b=>b.id===row.dataset.bodyResponseId),gene=S.config.genes?.[Number(row.dataset.bodyResponseGene)],preview=row.querySelector(".body-response-preview");if(body&&gene&&preview)preview.outerHTML=window.ChromatinModelExplain.bodyResponseHTML(body,gene.body_responses?.[body.id]);}); }
  function markDirty({preservePreset=false}={}) { pause();updateBodyResponsePreviews();S.dirty=true;if(!preservePreset){S.presetId="";$("#preset-select").value="";}updateControls();updateParticleCount();genomeView?.setDirty(true);interactionView?.setDirty(true);regulationView?.setDirty(true); }
  function updateParticleCount() {
    if(!S.config) return;
    const n=(Number(S.config.chromatin.beads)||0)+S.config.proteins.reduce((sum,p)=>sum+(Number(p.count)||0)*(p.representation==="idr"?Number(p.idr_beads)||1:1),0);
    const limit=S.capabilities?.limits?.particles;
    $("#particle-budget").textContent=`${pretty(n,0)}${limit?` / ${pretty(limit,0)}`:""} particles`;
    $("#particle-budget").classList.toggle("over-limit",!!limit&&n>limit);
    updateGenomicDisplays();
    const wall=$("#confinement-explanation");if(wall)wall.innerHTML=window.ChromatinModelExplain.confinementHTML(S.config,S.snapshot,S.applied);
  }
  function updateGenomicDisplays(){
    if(!S.config)return;
    const draft=genomeScale(S.config),applied=genomeScale(S.applied);
    $("#draft-genomic-scale").classList.toggle("dirty",S.dirty);
    $("#draft-genomic-scale").innerHTML=`<span>${S.dirty?"DRAFT GENOME":"GENOMIC MODEL"}</span><strong>${genomicSize(draft.total)}</strong><small>${pretty(draft.beads,0)} beads × ${genomicSize(draft.bp)} per bead${S.dirty?" · unapplied":""}</small>`;
    if(S.applied)$("#applied-genomic-scale").textContent=`Applied model · ${genomicSize(applied.total)} · ${pretty(applied.beads,0)} beads · ${genomicSize(applied.bp)}/bead`;
    const resolution=$("#chromatin-resolution-summary");
    if(resolution){resolution.innerHTML=`<strong>${genomicSize(draft.total)}</strong> total modeled span<small>${pretty(draft.beads,0)} beads × ${genomicSize(draft.bp)} per bead</small>`;}
    $$("[data-gene-genomic]").forEach(el=>{const gene=S.config.genes[Number(el.dataset.geneGenomic)];if(gene)el.textContent=geneGenomicSummary(gene,draft.bp);});
    $$("[data-regulatory-genomic]").forEach(el=>{const region=S.config.regulatory_regions[Number(el.dataset.regulatoryGenomic)];if(region)el.textContent=`${genomicSize(region.start*draft.bp)}–${genomicSize((region.end+1)*draft.bp)} · ${pretty(region.end-region.start+1,0)} beads`;});
    const planner=$("#nucleosome-resolution-plan");if(planner){const limit=S.capabilities?.limits?.chromatin_beads||4000;planner.innerHTML=`<strong>${genomicSize(draft.total)}</strong> at ${genomicSize(draft.bp)} per bead<small>${pretty(draft.beads,0)} current beads. At 200 bp per bead, this span needs ${pretty(Math.ceil(draft.total/200),0)} beads.</small><small>20 Mb at 200 bp needs <b>100,000 beads</b>; the local limit is ${pretty(limit,0)}. The 80 kb nucleosome-scale preset uses 400 beads.</small><small>200 bp labels one coarse nucleosome unit. Wrapped DNA and histone cores are not resolved. Optional eviction/reassembly is a binary occupancy model.</small>`;}
    const speed=$("#cohesin-speed-summary");
    if(speed){const arm=Number(S.config.cohesin?.step_rate||0)*draft.bp/1000;speed.textContent=`Nominal speed: ${pretty(arm,2)} kb per arm / reduced time. Two-arm loop growth: ${pretty(arm*2,2)} kb / reduced time before stalling.`;}
  }
  function updateControls() {
    $("#run-btn").disabled=!S.ready||S.dirty||S.busy&&!S.running;
    $("#step-btn").disabled=!S.ready||S.dirty||S.busy||S.running;
    $("#reset-btn").disabled=!S.ready||S.busy;
    $("#apply-btn").disabled=!S.config||!S.dirty||S.busy;
    $("#export-btn").disabled=!S.config;
    $("#export-run-btn").disabled=!S.snapshot||S.busy;
    $("#cluster-btn").disabled=!S.config;
    $("#learn-btn").disabled=!S.config;
    $("#preset-select").disabled=S.busy;
    $("#engine-select").disabled=S.busy;
    $("#run-btn").classList.toggle("running",S.running);
    $("#run-btn").innerHTML=S.running?'<span aria-hidden="true">Ⅱ</span> Pause':'<span aria-hidden="true">▶</span> Run simulation';
    $("#run-status").classList.toggle("running",S.running);
    $("#run-status").innerHTML=`<i></i> ${S.running?"Running":S.busy?"Updating":"Paused"}`;
    $("#draft-status").textContent=S.dirty?"Unapplied changes":S.ready?"All changes applied":"Connecting…";
    $("#draft-status").classList.toggle("dirty",S.dirty);
    const mapping=S.config?.calibration;
    $("#calibration-status").textContent=mapping?`Draft clock: ${mapping.status.replace(/_/g," ")} · Inspect ↗`:"No physical clock · Calibrate ↗";
    genomeView?.setRunState({running:S.running,busy:S.busy,ready:S.ready,dirty:S.dirty});
  }
  function pause() { S.running=false;runGeneration++;updateControls(); }
  async function settleRun() { pause();if(S.inflight)await S.inflight.catch(()=>{}); }
  async function reset(config=S.config, options={}) {
    const presetIdForReset=config===S.applied?(S.appliedPresetId||""):S.presetId;
    await settleRun(); S.busy=true;updateControls();$("#error-banner").hidden=true;
    try {
      const validation=await api("/api/validate",config);
      if(!validation.valid) throw new Error((validation.errors||["Invalid configuration"]).join("\n"));
      const snapshot=await api("/api/reset",{config,backend:S.backend});
      S.validationWarnings=validation.warnings||[];S.config=normalizeConfig(clone(config));S.applied=clone(S.config);S.appliedBackend=S.backend;S.presetId=presetIdForReset;S.appliedPresetId=presetIdForReset;$("#preset-select").value=S.presetId;S.dirty=false;S.ready=true;S.history=[];S.events=[];S.eventKeys.clear();S.selected=null;
      if(options.resetCamera) S.camera={yaw:-0.42,pitch:0.22,zoom:1};
      acceptSnapshot(snapshot);renderEditor();renderLegend();updateParticleCount();
      $("#connection-status").innerHTML='<i></i> Local engine connected';
      $("#experiment-title").textContent=S.applied.name||"Untitled experiment";
    } finally {S.busy=false;updateControls();}
  }
  async function stepOnce() {
    if(S.busy||!S.ready||S.dirty) return;
    S.busy=true;updateControls();
    try { S.inflight=api("/api/step",{steps:Number(S.applied.simulation.steps_per_frame)||10});acceptSnapshot(await S.inflight); }
    finally {S.inflight=null;S.busy=false;updateControls();}
  }
  async function run() {
    if(S.running) {pause();return;}
    if(S.dirty||!S.ready||S.busy) return;
    S.running=true;const generation=++runGeneration;updateControls();
    while(S.running && generation===runGeneration) {
      try {await stepOnce();} catch(e){pause();showError(e.message);break;}
      await new Promise(resolve=>setTimeout(resolve,45));
    }
  }
  function acceptSnapshot(snapshot) {
    S.snapshot=snapshot;
    const m=snapshot.metrics||{}, genes=snapshot.genes||[];
    const rna=genes.reduce((a,g)=>{a.nascent+=Number(g.nascent)||0;a.unprocessed+=Number(g.unprocessed)||0;a.mature+=Number(g.mature)||0;return a;},{nascent:0,unprocessed:0,mature:0});
    if(!S.history.length||S.history.at(-1).step!==snapshot.step) S.history.push({time:snapshot.time,step:snapshot.step,...m,...rna});
    if(S.history.length>1000) S.history=S.history.filter((_,i)=>i%2===0||i>S.history.length-200);
    for(const e of snapshot.events||[]) {
      const key=`${e.time}|${e.type}|${e.detail}`;
      if(!S.eventKeys.has(key)){S.eventKeys.add(key);S.events.push(e);}
    }
    S.events=S.events.slice(-160);
    $("#scene-empty").hidden=true;
    $("#metric-time").textContent=pretty(snapshot.time,2);
    $("#metric-step").textContent=`${pretty(snapshot.step,0)} integration steps`;
    const mapping=S.applied?.calibration;
    $("#metric-physical-time").textContent=mapping?.status==="fitted"&&Number.isFinite(mapping.time_s)?`${pretty(snapshot.time*mapping.time_s,3)} mapped seconds`:"No accepted physical clock";
    $("#metric-bound").textContent=`${pretty((m.bound_fraction||0)*100,1)}%`;
    $("#metric-cluster").textContent=pretty(m.largest_cluster,0);
    $("#metric-marked").textContent=`${pretty((m.marked_fraction||0)*100,1)}%`;
    $("#metric-rg").textContent=pretty(m.rg,2);
    $("#sample-count").textContent=`${S.history.length} samples`;
    $("#fraction-chart-value").textContent=`${pretty((m.bound_fraction||0)*100,0)}% bound`;
    $("#structure-chart-value").textContent=`Rg ${pretty(m.rg,2)}`;
    $("#rna-chart-value").textContent=`${pretty(rna.mature,0)} processed`;
    $("#engine-badge").textContent=snapshot.engine||S.backend;
    const warnings=[...new Set([...(S.validationWarnings||[]),...(snapshot.warnings||[])])];
    $("#runtime-warnings").hidden=!warnings.length;
    $("#runtime-warnings").innerHTML=warnings.length?`<details><summary>${warnings.length} model notes · reduced units and interpretation</summary><ul>${warnings.map(w=>`<li>${esc(w)}</li>`).join("")}</ul></details>`:"";
    genomeView?.update(snapshot,S.applied,S.selected,{dirty:S.dirty,running:S.running,busy:S.busy,ready:S.ready});
    interactionView?.update(snapshot,S.applied,{dirty:S.dirty});
    regulationView?.update(snapshot,S.applied,S.selected,{dirty:S.dirty});
    renderGenes();renderEvents();renderInspector();const wall=$("#confinement-explanation");if(wall)wall.innerHTML=window.ChromatinModelExplain.confinementHTML(S.config,S.snapshot,S.applied);S.chartsDirty=true;queueRender();
  }

  // Editors write only to the draft. Applying validates and restarts the engine.
  function field(path,label,opts={}) {
    const fieldNotes={"Free diffusion D":"Free-molecule mobility in length² / reduced time. Bound motion also depends on mechanical links.","Binding valence":"Maximum simultaneous protein and chromatin links per molecule.","Persistence length":"Resistance to bending, in reduced length units. A larger value favors straighter segments.","Integration time step":"Time advanced per step. Check smaller values before interpreting a result.","Random seed":"Sets stochastic choices. Independent replicas use different seeds."};
    opts={...opts,note:opts.note??fieldNotes[label]};
    const stored=getPath(path),value=stored===undefined?opts.defaultValue:stored,type=opts.type||"number",id=`f-${path.join("-").replace(/[^\w-]/g,"_")}`;
    const title=opts.help?` title="${esc(opts.help)}"`:"";
    const attr=`data-path="${pathAttr(path)}" id="${id}"${title}`;
    let control;
    if(type==="select") control=`<select ${attr}>${(opts.options||[]).map(o=>{const [v,n]=Array.isArray(o)?o:[o,o];return `<option value="${esc(v)}" ${String(value??"")===String(v)?"selected":""}>${esc(n)}</option>`;}).join("")}</select>`;
    else if(type==="checkbox") control=`<input ${attr} type="checkbox" ${value?"checked":""}>`;
    else control=`<input ${attr} type="${type}" value="${esc(value??(type==="number"?0:""))}" ${type==="number"?`step="${opts.step??"any"}" ${opts.min!==undefined?`min="${opts.min}"`:""} ${opts.max!==undefined?`max="${opts.max}"`:""}`:""}>`;
    return `<label for="${id}" class="form-field ${opts.full?"full":""} ${type==="checkbox"?"field-checkbox":""}"><span>${esc(label)}</span>${control}${opts.note?`<small>${esc(opts.note)}</small>`:""}</label>`;
  }
  const num=(path,label,opts={})=>field(path,label,{min:0,...opts});
  const textField=(path,label,opts={})=>field(path,label,{type:"text",...opts});
  function section(title,count="") {return `<div class="editor-subheading"><h3>${esc(title)}</h3><span>${esc(count)}</span></div>`;}
  function desc(title,description) {return `<h2>${esc(title)}</h2><p class="editor-description">${description}</p>`;}
  function addButton(action,label) {return `<button class="add-button" data-action="${action}">＋ ${esc(label)}</button>`;}
  function card(kind,index,name,color,meta,body,open=false) {return `<details class="entity-card" data-card="${kind}-${index}" ${open?"open":""}><summary><i class="entity-swatch" style="background:${esc(color||"#9aa98e")}"></i><span class="entity-name">${esc(name)}</span><span class="entity-meta">${esc(meta)}</span><button class="remove-button" data-action="remove-${kind}" data-index="${index}" aria-label="Remove ${esc(name)}" title="Remove">×</button></summary><div class="entity-body">${body}</div></details>`;}
  function affinityFields(path,entries,signed=false) {return `<div class="affinity-grid">${entries.map(([id,name,color])=>{
    const current=getPath([...path,id]),inherited=path.includes("chromatin_affinity")&&id!=="unmodified"&&current===undefined;
    return `<label class="affinity-label" for="aff-${path.join("-")}-${esc(id)}"><i class="affinity-dot" style="background:${esc(color||"#9ca7a1")}"></i>${esc(name)}</label><input id="aff-${path.join("-")}-${esc(id)}" aria-label="Affinity for ${esc(name)}" title="${signed?"Signed imposed-field affinity: negative repels, positive attracts.":inherited?"Unset: no override for this mark. Zero overrides the baseline; another present mapped mark can still permit binding.":"Dimensionless encounter-rate multiplier; clear a mark value to inherit the unmodified baseline."}" type="number" min="${signed?-20:0}" max="20" step="0.1" value="${esc(current??(inherited?"":0))}" placeholder="${inherited?"baseline":""}" data-path="${pathAttr([...path,id])}">`;
  }).join("")}</div>`;}

  function renderEditor() {
    if(["cluster","learn","calibration","studies"].includes(S.tab))S.builderExpanded=false;
    updateBuilderLayout();
    if(!S.config)return;
    const leavingCluster=(document.body.classList.contains("cluster-workspace")||document.body.classList.contains("learn-workspace")||document.body.classList.contains("calibration-workspace")||document.body.classList.contains("study-workspace"))&&!["cluster","learn","calibration","studies"].includes(S.tab);
    document.body.classList.toggle("cluster-workspace",S.tab==="cluster");
    document.body.classList.toggle("learn-workspace",S.tab==="learn");
    document.body.classList.toggle("calibration-workspace",S.tab==="calibration");
    document.body.classList.toggle("study-workspace",S.tab==="studies");
    $$("#editor-nav [data-tab]").forEach(b=>{b.setAttribute("aria-current",b.dataset.tab===S.tab?"page":"false");});
    if(clusterView){clusterView.destroy();clusterView=null;}
    if(learnView){learnView.destroy();learnView=null;}
    if(calibrationView){calibrationView.destroy();calibrationView=null;}
    if(studyView){studyView.destroy();studyView=null;}
    if(S.tab==="studies"){S.studySession??={};studyView=window.ChromatinStudy.render($("#editor-content"),{getConfig:()=>S.config,session:S.studySession,onNavigate:section=>navigateEditor(section)});return;}
    if(S.tab==="calibration"){
      S.calibrationSession??={};
      calibrationView=window.ChromatinCalibration.render($("#editor-content"),{config:S.config,profile:S.config.calibration,session:S.calibrationSession,onApply:profile=>{if(profile)S.config.calibration=clone(profile);else delete S.config.calibration;markDirty();toast(profile?"Physical mapping attached to the draft. Export it or apply to a fresh run.":"Mapping removed from the draft.");},onNavigate:section=>navigateEditor(section)});
      return;
    }
    if(S.tab==="learn"){
      learnView=window.ChromatinLearn.render($("#editor-content"),{config:S.config,initialTopic:S.learnTopic||"overview",onNavigate:section=>navigateEditor(section)});
      return;
    }
    if(leavingCluster)queueRender();
    if(S.tab==="cluster"){
      S.clusterJob??={...window.ChromatinClusterHelpers.DEFAULT_JOB};
      clusterView=new window.ChromatinClusterView($("#editor-content"),{job:S.clusterJob,getConfig:()=>S.config,onDownload:()=>toast("Cluster ZIP exported. Your local run is unchanged."),onExit:()=>{$("#editor-nav [data-tab='system']").click();$("#mobile-workspace-nav [data-mobile-tab='scene']").click();}});
      return;
    }
    const c=S.config;
    const openCards=new Set($$("#editor-content details[open]").map(el=>el.dataset.card));
    const scroll=$("#editor-content").scrollTop;
    let html="";
    if(S.tab==="system") {
      const preset=S.presets.find(p=>p.id===S.presetId);
      html=desc("The experiment","Define the simulation volume and integration settings. Inputs and trajectories use <strong>reduced units</strong>. Open <strong>Calibrate</strong> to attach a documented physical length and time mapping.");
      if(preset)html+=`<div class="preset-description">${esc(preset.description)}</div>`;
      html+=regulatoryModelSummary(c);
      html+=`<div class="form-grid">${textField(["name"],"Experiment name",{full:true})}${num(["seed"],"Random seed",{step:1})}${num(["simulation","dt"],"Integration time step",{min:0.000001,step:0.001})}${num(["simulation","steps_per_frame"],"Steps per display frame",{min:1,step:1,help:"Each sequential engine request advances this many integration steps."})}</div>`;
      html+=section("Nuclear confinement","soft wall");
      html+=`<div class="form-grid">${num(["simulation","confinement_radius"],"Reference nuclear radius R",{min:1,note:"Reduced length; initial/reference sphere, not a rigid membrane."})}${num(["simulation","confinement_stiffness"],"Wall stiffness k",{min:0,max:200,defaultValue:30,note:"0 disables the wall force. Legacy models use 30."})}</div><div id="confinement-explanation">${window.ChromatinModelExplain.confinementHTML(c,S.snapshot,S.applied)}</div>`;
      html+=`<div class="info-note">Smaller time steps improve stability. Display-frame size changes how often observables are sampled; it does not set physical diffusion.</div>`;
      html+=section("What is coupled");
      html+=`<ul class="model-list"><li>Chromatin motion and local stiffness</li><li>Diffusing proteins and coarse IDR chains</li><li>Finite-valence reversible binding</li><li>Occupancy-dependent chromatin modification</li><li>Chromatin-dependent gene activity</li><li>Imposed speckle and nucleolar fields</li></ul>`;
      html+=`<div class="info-note purple">Edit any section, then <strong>Apply & restart</strong>. Changes are validated before a new trajectory is created.</div>`;
    }
    if(S.tab==="proteins") {
      html=desc("Protein species","Mix freely diffusing particles with explicit coarse domain-bead chains. Binding changes mobility through mechanical links.");
      c.proteins.forEach((p,i)=>{
        const path=["proteins",i];
        let body=`<div class="form-grid">${textField([...path,"name"],"Species name",{full:true})}${field([...path,"color"],"Color",{type:"color"})}${num([...path,"count"],"Molecules",{step:1})}${num([...path,"diffusion"],"Free diffusion D")}${num([...path,"radius"],"Particle radius",{min:0.01})}${num([...path,"valence"],"Binding valence",{step:1})}${field([...path,"representation"],"Representation",{type:"select",options:[["bead","Single bead"],["idr","Coarse IDR chain"]]})}</div>`;
        const roles=proteinRoles(c,p.id);if(roles.length)body=`<div class="info-note"><strong>Configured role:</strong> ${esc(roles.join(" · "))}. This is an explicit diffusing species; its role follows its assignments in Genes and Regulation.</div>`+body;
        if(p.representation==="idr")body+=section("Coarse domain chain")+`<div class="form-grid">${num([...path,"idr_beads"],"Beads per molecule",{step:1,min:2})}${num([...path,"sticker_fraction"],"Sticker fraction",{max:1,step:0.05})}</div><div class="info-note purple">A domain-level chain with stickers, not a residue-level CALVADOS force field. Valence is shared by chromatin and protein links.</div>`;
        body+=section("Chromatin encounter affinity")+affinityFields([...path,"chromatin_affinity"],[["unmodified","Unmodified","#869cb0"],...c.marks.map(m=>[m.id,m.name,m.color])]);
        body+=`<p class="editor-description" style="font-size:9px;margin:9px 0 0">Dimensionless association-rate multipliers. “Baseline” means no override. Present mapped marks use their greatest value, including explicit 0; otherwise the unmodified baseline applies. Another positive mapped mark can still permit binding.</p>`;
        body+=section("Nucleosome access & release")+`<div class="form-grid">${num([...path,"nucleosome_affinity","occupied"],"Occupied-site capture ×",{note:"Association multiplier when this site has a nucleosome; 0 blocks new binding."})}${num([...path,"nucleosome_affinity","vacant"],"Vacant-site capture ×",{note:"Association multiplier in a nucleosome-free site; species may prefer either state."})}${num([...path,"nucleosome_off_rate_multiplier","occupied"],"Occupied-site release ×",{note:"Multiplies the unbinding hazard while the bound site is occupied."})}${num([...path,"nucleosome_off_rate_multiplier","vacant"],"Vacant-site release ×",{note:"Multiplies the unbinding hazard while the bound site is vacant; 1 preserves the base rate."})}</div><div class="info-note">Capture hazard = binding rate × mark preference × element preference × occupancy preference, conditional on an eligible encounter. Release hazard = unbinding rate × current occupancy release factor. These are reduced-time hazards, not molar association constants.</div>`;
        body+=section("Chromatin remodeling","while recruited")+`<div class="form-grid">${num([...path,"nucleosome_eviction_rate"],"Nucleosome eviction rate",{note:"Remove a nucleosome and its marks at a recruited site."})}${num([...path,"nucleosome_reassembly_rate"],"Nucleosome reassembly rate",{note:"Fill a vacant recruited site with an unmarked nucleosome."})}${num([...path,"nucleosome_sliding_rate"],"Nucleosome sliding rate",{full:true,note:"Move an occupied nucleosome and its marks into an adjacent vacant bead."})}</div><p class="editor-description" style="font-size:10px;margin:9px 0 0">Requires enabled nucleosome occupancy and chromatin recruitment. Set a rate to 0 to disable that action. Sliding preserves nucleosome number. These coarse kinetic rules do not resolve ATP hydrolysis or DNA wrapping.</p>`;
        body+=section("Binding & catalysis")+`<div class="form-grid">${num([...path,"on_rate"],"Binding rate")}${num([...path,"off_rate"],"Unbinding rate")}${field([...path,"writer_mark"],"Write mark",{type:"select",options:[["","None"],...c.marks.map(m=>[m.id,m.name])]})}${field([...path,"eraser_mark"],"Erase mark",{type:"select",options:[["","None"],...c.marks.map(m=>[m.id,m.name])]})}${num([...path,"catalytic_rate"],"Catalytic rate",{full:true,note:"Modification requires binding to chromatin."})}</div>`;
        html+=card("protein",i,p.name,p.color,`${p.count} × ${p.representation==="idr"?"chain":"bead"}`,body,i===0);
      });
      html+=addButton("add-protein","Add protein species");
      html+=addButton("add-remodeler","Add chromatin remodeler");
    }
    if(S.tab==="interactions") {
      html=desc("Protein interactions","Specify reversible links between species. Strength sets link forces; association and dissociation rates set their kinetics.");
      if(!c.proteins.length)html+='<div class="info-note">Add protein species first to configure their interactions.</div>';
      else {
        html+=section("Link strength matrix","symmetric");
        html+=`<div class="matrix-wrap"><table class="interaction-matrix"><thead><tr><th></th>${c.proteins.map(p=>`<th title="${esc(p.name)}">${esc(p.name.slice(0,6))}</th>`).join("")}</tr></thead><tbody>${c.proteins.map((p,i)=>`<tr><th class="row-head" title="${esc(p.name)}">${esc(p.name.slice(0,7))}</th>${c.proteins.map((q,j)=>{const pair=c.interactions.find(x=>(x.a===p.id&&x.b===q.id)||(x.a===q.id&&x.b===p.id));return `<td><input type="number" min="0" step="0.1" aria-label="${esc(p.name)} to ${esc(q.name)} interaction strength" data-pair-a="${esc(p.id)}" data-pair-b="${esc(q.id)}" value="${pair?.strength||0}" class="${pair?.strength?"has-interaction":""} ${j<i?"mirrored":""}"></td>`;}).join("")}</tr>`).join("")}</tbody></table></div>`;
        html+=`<div class="info-note">Nonzero entries create link rules. Setting spring strength to 0 does not stop association: set the association rate to 0 or remove the rule. Bonds share each molecule’s total valence. Inspect live connections in the <strong>Interactions</strong> simulation view.</div>`;
        html+=section("Pair kinetics",`${c.interactions.length} rules`);
        c.interactions.forEach((p,i)=>{const a=c.proteins.find(x=>x.id===p.a),b=c.proteins.find(x=>x.id===p.b),path=["interactions",i];html+=card("interaction",i,`${a?.name||p.a} ↔ ${b?.name||p.b}`,a?.color,"",`<div class="form-grid">${num([...path,"strength"],"Link strength")}${num([...path,"radius"],"Encounter radius",{min:.01})}${num([...path,"on_rate"],"Association rate")}${num([...path,"off_rate"],"Dissociation rate")}</div>`,i===0);});
      }
    }
    if(S.tab==="chromatin") {
      html=desc("Chromatin polymer","Set genomic span, resolution, polymer mechanics, and independent marks. The 20 Mb preset uses 2,000 beads at 10 kb per bead.");
      html+='<div id="chromatin-resolution-summary" class="genome-resolution-summary" aria-live="polite"></div>';
      html+='<details class="resolution-planner"><summary>Nucleosome-scale resolution planner</summary><div id="nucleosome-resolution-plan" aria-live="polite"></div></details>';
      html+=`<div class="form-grid">${field(["chromatin","representation"],"Chromatin representation",{type:"select",full:true,options:[["coarse","Coarse genomic beads"],["nucleosome_unit","Coarse nucleosome units"]]})}${num(["chromatin","beads"],"Chromatin beads",{min:8,max:S.capabilities?.limits?.chromatin_beads||4000,step:1})}${num(["chromatin","bp_per_bead"],"Base pairs per bead",{min:1,max:1000000,step:1})}${num(["chromatin","bond_length"],"Bond length",{min:.01})}${num(["chromatin","bond_k"],"Bond stiffness",{min:.01})}${num(["chromatin","diffusion"],"Bead diffusion D")}${num(["chromatin","persistence_length"],"Persistence length")}${num(["chromatin","repulsion"],"Excluded-volume strength",{full:true})}</div>`;
      html+=`<div class="info-note">Local limits: <strong>${pretty(S.capabilities?.limits?.chromatin_beads||4000,0)} chromatin beads</strong> and <strong>${pretty(S.capabilities?.limits?.particles||6000,0)} total particles</strong>. Span = bead count × base pairs per bead.</div><div class="info-note purple">Manual resolution changes keep gene, regulatory-element, mark-region, and barrier bead indices unchanged. Review their genomic positions and endpoints after editing. Mechanics and reaction rates are not automatically rescaled.</div>`;
      html+=section("Nucleosome occupancy","optional");
      const nucleosomePath=["chromatin","nucleosomes"];
      html+=`<div class="form-grid">${field([...nucleosomePath,"enabled"],"Enable nucleosome occupancy",{type:"checkbox"})}${num([...nucleosomePath,"initial_occupancy"],"Initial occupied fraction",{max:1,step:.05,note:"Initial occupancy probability; intervals below override it."})}${num([...nucleosomePath,"eviction_rate"],"Background eviction rate")}${num([...nucleosomePath,"reassembly_rate"],"Background reassembly rate")}${num([...nucleosomePath,"occupied_persistence_multiplier"],"Occupied persistence ×")}${num([...nucleosomePath,"vacant_persistence_multiplier"],"Vacant persistence ×")}</div>`;
      html+='<div class="info-note purple">One binary state per chromatin bead. At coarse genomic resolution this is an effective occupancy flag, not one physical nucleosome. Eviction clears marks; vacant beads cannot acquire histone marks. Reassembly starts unmarked. Persistence factors combine with transcription-dependent stiffness.</div>';
      html+=section("Initial occupied / vacant intervals");
      c.chromatin.nucleosomes.initial_regions.forEach((r,i)=>{const path=[...nucleosomePath,"initial_regions",i];html+=`<div class="range-title">Occupancy interval ${i+1}<button class="remove-button" data-action="remove-nucleosome-region" data-index="${i}" aria-label="Remove nucleosome occupancy interval ${i+1}">×</button></div><div class="form-grid">${num([...path,"start"],"Start bead",{step:1})}${num([...path,"end"],"End bead",{step:1})}${field([...path,"occupied"],"Initially occupied",{type:"checkbox"})}</div>`;});
      html+=addButton("add-nucleosome-region","Add nucleosome-free interval");
      html+='<p class="editor-description" style="font-size:10px">New intervals start vacant. These set the initial state; reassembly and remodeling can change it during a run. Occupancy changes each species’ capture and optional release hazards immediately at the next reaction update. Inspect individual sites in <strong>1D genome → Live binding</strong>.</p>';
      html+=section("State-selective chromatin attraction",`${c.chromatin.state_interactions.length} rules`);
      const chromatinStates=[["unmodified","Unmodified"],...c.marks.map(m=>[m.id,m.name])];
      c.chromatin.state_interactions.forEach((rule,i)=>{const path=["chromatin","state_interactions",i],a=chromatinStates.find(state=>state[0]===rule.a)?.[1]||rule.a,b=chromatinStates.find(state=>state[0]===rule.b)?.[1]||rule.b;html+=card("state-interaction",i,`${a} ↔ ${b}`,"#8c70b4","",`<div class="form-grid">${field([...path,"a"],"Chromatin state A",{type:"select",options:chromatinStates})}${field([...path,"b"],"Chromatin state B",{type:"select",options:chromatinStates})}${num([...path,"strength"],"Well depth · kBT",{max:20})}${num([...path,"radius"],"Attraction radius",{min:.2,max:5})}${num([...path,"min_separation"],"Minimum contour separation",{min:1,max:4000,step:1,full:true,note:"Exclude bead pairs closer along the polymer contour."})}</div>`,i===0);});
      html+=addButton("add-state-interaction","Add chromatin attraction rule");
      html+='<div class="info-note purple">These rules attract chromatin beads with the chosen marks, independently of protein binding. At each eligible bead pair, only the well with the most negative current energy applies. Live contacts and energies appear in Regulation.</div>';
      html+=section("Independent mark types",`${c.marks.length} marks`);
      c.marks.forEach((m,i)=>{const path=["marks",i];html+=card("mark",i,m.name,m.color,"",`<div class="form-grid">${textField([...path,"name"],"Mark name",{full:true})}${field([...path,"color"],"Color",{type:"color"})}${num([...path,"decay_rate"],"Decay rate")}${num([...path,"spread_rate"],"Contact-spread rate")}${num([...path,"contact_radius"],"Contact radius",{min:.01})}</div>`,i===0);});
      html+=addButton("add-mark","Add chromatin mark");
      html+='<div class="info-note purple">Marks can coexist on a bead. Contact-mediated spreading is an optional phenomenological rule, separate from explicit protein writers.</div>';
      html+=section("Initial marked regions");
      c.chromatin.initial_marks.forEach((r,i)=>{const path=["chromatin","initial_marks",i];html+=`<div class="range-title">Region ${i+1}</div><div class="range-controls">${num([...path,"start"],"Start",{step:1})}${num([...path,"end"],"End",{step:1})}${field([...path,"mark"],"Mark",{type:"select",options:c.marks.map(m=>[m.id,m.name])})}<button class="remove-button" data-action="remove-region" data-index="${i}" aria-label="Remove initial mark region ${i+1}">×</button></div>`;});
      html+=addButton("add-region","Add marked region");
    }
    if(S.tab==="genes") {
      html=desc("Transcription & RNA","Stochastic promoter ON/OFF switching and initiation respond to chromatin and bound transcription factors. Elongation, RNA processing and degradation are separate processes.");
      html+=regulatoryModelSummary(c);
      html+=section("Individual RNA molecules","optional");
      html+=`<div class="form-grid">${field(["transcription","enabled"],"Track individual transcripts",{type:"checkbox"})}${num(["transcription","max_transcripts"],"Maximum tracked transcripts",{min:1,max:5000,step:1,full:true})}${num(["transcription","rna_diffusion"],"Released RNA diffusion D")}${num(["transcription","mature_rna_diffusion"],"Processed RNA diffusion D")}</div><div class="info-note">Each RNA is one molecule with a growing nucleotide length, then a diffusing point after release. RNA chains and RNA–protein interactions are not resolved. The Regulation view shows polymerase progression, pauses, and transcript state. The RNA molecule cap is separate from the chromatin/protein particle budget.</div>`;
      c.genes.forEach((g,i)=>{const path=["genes",i];let body=`<div class="form-grid">${textField([...path,"name"],"Locus name",{full:true})}${num([...path,"start"],"Start bead",{step:1})}${num([...path,"end"],"End bead",{step:1})}</div>`;
        body+=`<div class="gene-genomic-summary" data-gene-genomic="${i}">${esc(geneGenomicSummary(g,Number(c.chromatin.bp_per_bead)))}</div>`;
        body+=section("Promoter & polymerase kinetics")+`<div class="form-grid">${num([...path,"on_rate"],"Promoter ON rate")}${num([...path,"off_rate"],"Promoter OFF rate")}${num([...path,"initiation_rate"],"Initiation rate")}${num([...path,"elongation_rate"],"Elongation · beads / time")}${num([...path,"processing_rate"],"RNA processing rate")}${num([...path,"rna_decay_rate"],"RNA decay rate")}</div>`;
        body+=section("Elongation & nucleosome coupling")+`<div class="form-grid">${num([...path,"pause_rate"],"Polymerase pause rate")}${num([...path,"resume_rate"],"Polymerase resume rate")}${num([...path,"occupied_elongation_multiplier"],"Occupied elongation ×",{note:"Speed factor while the current bead is occupied."})}${num([...path,"nucleosome_eviction_rate"],"Polymerase eviction rate")}</div><p class="editor-description" style="font-size:9px;margin:9px 0 0">Pause/resume tracking uses individual transcript mode. Nucleosome-dependent speed and eviction require occupancy to be enabled.</p>`;
        body+=section("Local regulation")+`<div class="form-grid">${field([...path,"activator_species"],"Gene-wide transcription factor",{type:"select",full:true,note:"Used only when this gene has no enabled regulatory elements.",options:[["","None"],...c.proteins.map(p=>[p.id,p.name])]})}${field([...path,"required_mark"],"Activating mark",{type:"select",options:[["","None"],...c.marks.map(m=>[m.id,m.name])]})}${field([...path,"repressive_mark"],"Repressive mark",{type:"select",options:[["","None"],...c.marks.map(m=>[m.id,m.name])]})}${num([...path,"speckle_gain"],"Speckle processing gain",{note:"Enhances processing, not initiation."})}${num([...path,"persistence_multiplier"],"Active stiffness multiplier",{min:.01,note:"1 disables polymerase-dependent bending changes."})}</div>`;
        body+=`<details class="body-response-editor"><summary>Nuclear environment → promoter rates</summary><p>Optional causal hypotheses for this gene. A multiplier of 1 is neutral; below 1 represses and above 1 enhances. These effects are separate from speckle-dependent RNA processing.</p><div class="mechanism-equation">w = exp[−d² / (2R²)]<br>rate factor = exp[w · ln(m)]</div><p>Distance d is from the gene-start promoter bead to the field center. Responses from enabled fields multiply. Manuscript proximity/SIF associations do not fit or prove this rate law.</p>${c.bodies.map(b=>{const response=g.body_responses?.[b.id]||{};return `<div class="body-response-row" data-body-response-gene="${i}" data-body-response-id="${esc(b.id)}"><h4>${esc(b.name)} <small>${b.enabled===false?"disabled field":b.kind}</small></h4><div class="form-grid">${num([...path,"body_responses",b.id,"activation_multiplier"],"Activation at field center ×",{min:.01,max:20,defaultValue:1})}${num([...path,"body_responses",b.id,"initiation_multiplier"],"Initiation at field center ×",{min:.01,max:20,defaultValue:1})}</div>${window.ChromatinModelExplain.bodyResponseHTML(b,response)}</div>`;}).join("")||"<p>Add a nuclear field in Bodies to configure a response.</p>"}</details>`;
        html+=card("gene",i,g.name,"#59a68e",`${g.start}–${g.end}`,body,i===0);
      });
      html+=addButton("add-gene","Add gene locus");
      html+='<div class="info-note">Bead endpoints are inclusive. Genomic summary ranges show the outer boundaries of the covered interval. Regulatory couplings are explicit hypotheses; processed RNA, nascent polymerases, and unprocessed transcripts are distinct observables.</div>';
    }
    if(S.tab==="regulation") {
      html=desc("Enhancers & promoters","Place regulatory elements along chromatin, choose their transcription factor and target gene, then define contact-dependent transcriptional gains.");
      html+='<div class="info-note">TF affinity boosts act only inside enabled elements. Coactivators are recruited by actual protein–protein links to DNA-bound TFs; configure those link rules in Interactions and mark-writing activity in Proteins.</div>';
      c.regulatory_regions.forEach((r,i)=>{const path=["regulatory_regions",i],body=`<div class="form-grid">${field([...path,"enabled"],"Enable element",{type:"checkbox"})}${textField([...path,"name"],"Element name",{full:true})}${field([...path,"kind"],"Element type",{type:"select",options:[["enhancer","Enhancer"],["promoter","Promoter"]]})}${field([...path,"color"],"Track color",{type:"color"})}${num([...path,"start"],"Start bead",{step:1})}${num([...path,"end"],"End bead",{step:1})}</div><div class="gene-genomic-summary" data-regulatory-genomic="${i}">${genomicSize(r.start*c.chromatin.bp_per_bead)}–${genomicSize((r.end+1)*c.chromatin.bp_per_bead)}</div>${section("Binding & target")}<div class="form-grid">${field([...path,"target_gene"],"Target gene",{type:"select",full:true,options:c.genes.map(g=>[g.id,g.name])})}${field([...path,"tf_species"],"Transcription factor",{type:"select",full:true,options:c.proteins.map(p=>[p.id,p.name])})}${field([...path,"coactivator_species"],"Recruited coactivator",{type:"select",full:true,options:[["","None"],...c.proteins.map(p=>[p.id,p.name])]})}${num([...path,"tf_affinity_multiplier"],"TF region affinity ×",{note:"Multiplies mark-dependent affinity; 1 adds no region preference."})}${num([...path,"contact_radius"],"Enhancer contact radius",{min:.01,note:"Reduced length; gates enhancer action on its target promoter."})}</div>${section("Transcription coupling")}<div class="form-grid">${num([...path,"activation_gain"],"Promoter activation gain")}${num([...path,"initiation_gain"],"Initiation gain")}${num([...path,"coactivator_gain"],"Coactivator signal gain",{full:true})}</div><div class="info-note purple">Gains are explicit model hypotheses. The Regulation view shows current TF binding, recruited molecules, contact gates, effective rates, and RNA output.</div>`;html+=card("regulatory",i,r.name,r.color,r.kind,body,i===0);});
      html+=`<div class="button-row">${addButton("add-enhancer","Add enhancer")}${addButton("add-promoter","Add promoter")}</div>`;
      html+='<p class="editor-description" style="margin-top:15px">Intervals include both bead endpoints. One enabled promoter per gene must contain its start bead. Enhancers must contact that promoter, or the gene start if no promoter is defined. Multiple enhancers can target the same gene.</p>';
    }
    if(S.tab==="bodies") {
      html=desc("Nuclear landmarks","These fixed spatial affinity fields represent speckle and nucleolar environments. Their positions and radii are imposed.");
      c.bodies.forEach((b,i)=>{const path=["bodies",i];let body=`<div class="form-grid">${field([...path,"enabled"],"Enable landmark",{type:"checkbox"})}${textField([...path,"name"],"Landmark name",{full:true})}${field([...path,"kind"],"Kind",{type:"select",options:[["speckle","Speckle field"],["nucleolus","Nucleolar field"]]})}${field([...path,"color"],"Color",{type:"color"})}${num([...path,"radius"],"Radius",{min:.01})}${num([...path,"strength"],"Field strength")}${field([...path,"chromatin_affinity"],"Chromatin affinity",{full:true,help:"Negative values repel chromatin; positive values attract. The server validates the supported range."})}</div>`;
        body+='<div class="mini-body-title">Position in reduced length units</div>'+`<div class="body-center-grid">${field([...path,"center",0],"X")}${field([...path,"center",1],"Y")}${field([...path,"center",2],"Z")}</div>`;
        body+=section("Protein partitioning affinity")+affinityFields([...path,"protein_affinity"],c.proteins.map(p=>[p.id,p.name,p.color]),true);
        html+=card("body",i,b.name,b.color,b.kind,body,i===0);
      });
      html+=`<div class="button-row">${addButton("add-speckle","Speckle field")}${addButton("add-nucleolus","Nucleolar field")}</div>`;
      html+='<div class="info-note purple"><strong>Imposed landmarks ≠ emergent condensates.</strong> Protein links can produce assemblies independently. Nucleolar RNA processing is not modeled.</div>';
      html+=section("Whole nucleus & tethering","envelope");
      if(!c.nucleus){html+=`<p class="editor-description">Add a nuclear envelope around the modeled territory, then define genomic tethering elements (for example lamina-associated domains) that attach to it or to a landmark surface.</p><div class="button-row">${addButton("add-nucleus","Whole nucleus")}</div>`;}
      else{const n=c.nucleus;html+=`<div class="form-grid">${field(["nucleus","enabled"],"Enable envelope",{type:"checkbox"})}${num(["nucleus","envelope_stiffness"],"Envelope stiffness k",{min:0,max:200,defaultValue:30,note:"Soft ellipsoidal wall; 0 removes its force."})}</div><div class="mini-body-title">Semi-axes (reduced length)</div><div class="body-center-grid">${field(["nucleus","radii",0],"a")}${field(["nucleus","radii",1],"b")}${field(["nucleus","radii",2],"c")}</div><div class="mini-body-title">Envelope center relative to the territory origin</div><div class="body-center-grid">${field(["nucleus","center",0],"X")}${field(["nucleus","center",1],"Y")}${field(["nucleus","center",2],"Z")}</div><div class="info-note">The territory sphere (System → confinement) still bounds the modeled domain. Place the envelope so it intersects that sphere where the domain should touch the nuclear periphery. The origin must stay inside the envelope.</div>`;
        (n.tethers||[]).forEach((t,i)=>{const path=["nucleus","tethers",i];let body=`<div class="form-grid">${field([...path,"enabled"],"Enable element",{type:"checkbox"})}${textField([...path,"name"],"Element name",{full:true})}${field([...path,"color"],"Track color",{type:"color"})}${field([...path,"anchor"],"Anchor",{type:"select",options:[["envelope","Nuclear envelope"],["body","Landmark surface"]]})}${num([...path,"start"],"Start bead",{step:1})}${num([...path,"end"],"End bead",{step:1})}${t.anchor==="body"?field([...path,"body_id"],"Landmark",{type:"select",full:true,options:c.bodies.map(b=>[b.id,b.name])}):""}</div>${section("Attachment kinetics")}<div class="form-grid">${num([...path,"attach_rate"],"Attach rate",{note:"Hazard per bead within capture distance."})}${num([...path,"detach_rate"],"Detach rate")}${num([...path,"capture_distance"],"Capture distance",{min:.1})}${num([...path,"spring_k"],"Tether spring k")}${num([...path,"rest_length"],"Rest length")}${field([...path,"required_mark"],"Required mark",{type:"select",options:[["","None"],...c.marks.map(m=>[m.id,m.name])]})}</div>${section("Nucleosome occupancy dependence")}<div class="form-grid">${num([...path,"occupancy_attach_multiplier","occupied"],"Occupied attach ×",{defaultValue:1})}${num([...path,"occupancy_attach_multiplier","vacant"],"Vacant attach ×",{defaultValue:1})}${num([...path,"occupancy_detach_multiplier","occupied"],"Occupied detach ×",{defaultValue:1})}${num([...path,"occupancy_detach_multiplier","vacant"],"Vacant detach ×",{defaultValue:1})}</div><div class="info-note purple">Phenomenological anchor springs, not lamina proteins or pore complexes. Rates are reduced-time hazards and model inputs.</div>`;html+=card("tether",i,t.name,t.color,`${t.start}–${t.end} · ${t.anchor==="body"?"landmark":"envelope"}`,body,i===0);});
        html+=`<div class="button-row">${addButton("add-tether","Tethering element")}<button class="secondary-button" data-action="remove-nucleus">Remove nucleus</button></div>`;}
    }
    if(S.tab==="cohesin") {
      html=desc("Loop extrusion","Cohesin pairs load on adjacent beads, step outward, stall at barriers or other extruders, and unload.");
      html+=`<div class="form-grid">${field(["cohesin","enabled"],"Enable cohesin extrusion",{type:"checkbox"})}${num(["cohesin","count"],"Extruder count",{step:1})}${num(["cohesin","step_rate"],"Extrusion step rate")}${num(["cohesin","off_rate"],"Unloading rate")}${num(["cohesin","spring_k"],"Loop spring stiffness")}</div>`;
      html+='<div id="cohesin-speed-summary" class="cohesin-speed-summary" aria-live="polite"></div>';
      html+=section("Directional barriers",`${c.chromatin.barriers.length} sites`);
      c.chromatin.barriers.forEach((b,i)=>{const path=["chromatin","barriers",i];html+=`<div class="range-title">Barrier ${i+1}</div><div class="range-controls" style="grid-template-columns:1fr 1.3fr 18px">${num([...path,"bead"],"Bead index",{step:1})}${field([...path,"direction"],"Blocks motion",{type:"select",options:[["left","Leftward"],["right","Rightward"],["both","Both directions"]]})}<button class="remove-button" data-action="remove-barrier" data-index="${i}" aria-label="Remove barrier ${i+1}">×</button></div>`;});
      html+=addButton("add-barrier","Add extrusion barrier");
      html+='<div class="info-note">Extrusion kinetics use reduced time. Amber links in the scene show the current loop endpoints.</div>';
    }
    if(S.tab==="model") {
      html=desc("About this model","An inspectable research prototype for spatial coupling, with explicit unit mapping and reproducible reference benchmarks.");
      for(const b of S.capabilities?.backends||[])html+=`<div class="backend-card"><strong>${esc(b.id==="numpy"?"Brownian reference":b.id==="polychrom"?"Polychrom / OpenMM":b.id)}</strong><p>${esc(b.detail)}</p><span class="backend-status ${b.available?"":"unavailable"}">${b.available?"Available locally":"Not installed"}</span></div>`;
      html+=section("Interpretation");
      html+='<ul class="model-list"><li>kBT = 1. Lengths and time are reduced units.</li><li>Binding physically couples particles; free diffusion alone is not an observed bound diffusion coefficient.</li><li>IDRs are coarse domain chains, not a sequence-calibrated force field.</li><li>Speckles and nucleoli are fixed fields. Explicit protein clusters are a separate mechanism.</li><li>Local stiffness and transcription couplings are configurable hypotheses.</li><li>Population averages and single trajectories require independent seeds and convergence checks.</li></ul>';
      if(S.capabilities?.limitations?.length)html+=section("Engine limitations")+`<ul class="model-list">${S.capabilities.limitations.map(x=>`<li>${esc(x)}</li>`).join("")}</ul>`;
      html+='<div class="info-note purple">Use <strong>Export</strong> to save a configuration. <strong>Save run</strong> includes the applied model, state, and server-recorded history.</div>';
    }
    const explanations={system:["motion","How motion and time steps work"],proteins:["binding","How proteins move, bind and assemble"],interactions:["binding","From encounters to binding affinity"],chromatin:["chromatin","How chromatin access and stiffness work"],genes:["transcription","From promoter switching to a transcript"],regulation:["transcription","How regulatory inputs change transcription"],bodies:["architecture","How nuclear compartments shape motion"],cohesin:["architecture","How loop extrusion works"],model:["overview","Explore the illustrated model guide"]};
    const help=explanations[S.tab];
    $("#editor-content").innerHTML=(help?`<button class="context-lesson" data-open-lesson="${help[0]}"><span aria-hidden="true">◈</span><span>${help[1]}</span><span aria-hidden="true">↗</span></button>`:"")+html;
    if(openCards.size)$$("#editor-content details").forEach(d=>d.open=openCards.has(d.dataset.card));
    $("#editor-content").scrollTop=scroll;
    updateGenomicDisplays();
  }

  function uniqueId(prefix) {return `${prefix}_${Date.now().toString(36)}_${Math.random().toString(36).slice(2,5)}`;}
  function handleAction(action,index) {
    const c=S.config;if(!c)return;
    if(action==="add-protein"||action==="add-remodeler") {
      const remodeling=action==="add-remodeler",id=uniqueId(remodeling?"remodeler":"protein");c.proteins.push({id,name:remodeling?"Chromatin remodeler":`Protein ${c.proteins.length+1}`,color:colors[c.proteins.length%colors.length],count:12,diffusion:0.7,radius:0.3,valence:2,representation:"bead",idr_beads:5,sticker_fraction:0.6,chromatin_affinity:{unmodified:0.3,...Object.fromEntries(c.marks.map(m=>[m.id,remodeling?.3:0]))},nucleosome_affinity:{occupied:1,vacant:1},nucleosome_off_rate_multiplier:{occupied:1,vacant:1},nucleosome_eviction_rate:remodeling?.2:0,nucleosome_reassembly_rate:0,nucleosome_sliding_rate:remodeling?.5:0,on_rate:0.4,off_rate:0.15,writer_mark:"",eraser_mark:"",catalytic_rate:0.2});
      if(remodeling){c.chromatin.nucleosomes.enabled=true;toast("Remodeler added and nucleosome occupancy enabled in the draft. Review rates and recruitment, then Apply & restart.");}
    } else if(action==="remove-protein") {
      const [p]=c.proteins.splice(index,1);c.interactions=c.interactions.filter(x=>x.a!==p.id&&x.b!==p.id);c.bodies.forEach(b=>delete b.protein_affinity[p.id]);c.genes.forEach(g=>{if(g.activator_species===p.id)g.activator_species="";});c.regulatory_regions=c.regulatory_regions.filter(r=>r.tf_species!==p.id);c.regulatory_regions.forEach(r=>{if(r.coactivator_species===p.id)r.coactivator_species="";});
    } else if(action==="add-mark") {const id=uniqueId("mark");c.marks.push({id,name:`Mark ${c.marks.length+1}`,color:colors[(c.marks.length+2)%colors.length],decay_rate:0.01,spread_rate:0,contact_radius:1.2});
    } else if(action==="remove-mark") {
      const [m]=c.marks.splice(index,1);c.chromatin.initial_marks=c.chromatin.initial_marks.filter(r=>r.mark!==m.id);c.proteins.forEach(p=>{delete p.chromatin_affinity[m.id];if(p.writer_mark===m.id)p.writer_mark="";if(p.eraser_mark===m.id)p.eraser_mark="";});c.genes.forEach(g=>{if(g.required_mark===m.id)g.required_mark="";if(g.repressive_mark===m.id)g.repressive_mark="";});c.chromatin.state_interactions=c.chromatin.state_interactions.filter(r=>r.a!==m.id&&r.b!==m.id);
    } else if(action==="add-gene") {const start=Math.max(0,Math.floor(c.chromatin.beads*.35));c.genes.push({id:uniqueId("gene"),name:`Gene ${c.genes.length+1}`,start,end:Math.min(c.chromatin.beads-1,start+10),on_rate:0.05,off_rate:0.1,initiation_rate:0.6,elongation_rate:2,processing_rate:0.2,rna_decay_rate:0.01,activator_species:"",required_mark:"",repressive_mark:"",speckle_gain:2,persistence_multiplier:1,nucleosome_eviction_rate:0,pause_rate:0,resume_rate:1,occupied_elongation_multiplier:1});
    } else if(action==="remove-gene"){const [g]=c.genes.splice(index,1);c.regulatory_regions=c.regulatory_regions.filter(r=>r.target_gene!==g.id);}
    else if(action==="add-enhancer"||action==="add-promoter") {if(!c.genes.length||!c.proteins.length){toast("Add a gene and a protein species before defining a regulatory element.");return;}if(c.regulatory_regions.length>=128){toast("This model already has the supported 128 regulatory elements.");return;}const kind=action==="add-promoter"?"promoter":"enhancer",gene=kind==="promoter"?c.genes.find(g=>!c.regulatory_regions.some(r=>r.kind==="promoter"&&r.enabled&&r.target_gene===g.id)):c.genes[0];if(!gene){toast("Every gene already has an enabled promoter. Edit or disable an existing promoter before adding another.");return;}const start=kind==="promoter"?gene.start:Math.max(0,gene.start-20);c.regulatory_regions.push({id:uniqueId(kind),name:`${kind==="promoter"?"Promoter":"Enhancer"} ${c.regulatory_regions.filter(r=>r.kind===kind).length+1}`,color:kind==="promoter"?"#639f9c":"#c7975c",kind,start,end:Math.min(c.chromatin.beads-1,start+2),target_gene:gene.id,tf_species:c.proteins[0].id,coactivator_species:"",tf_affinity_multiplier:3,activation_gain:2,initiation_gain:1,coactivator_gain:1,contact_radius:2,enabled:true});}
    else if(action==="remove-regulatory")c.regulatory_regions.splice(index,1);
    else if(action==="add-nucleosome-region")c.chromatin.nucleosomes.initial_regions.push({start:0,end:Math.min(4,c.chromatin.beads-1),occupied:false});
    else if(action==="remove-nucleosome-region")c.chromatin.nucleosomes.initial_regions.splice(index,1);
    else if(action==="add-state-interaction"){const states=["unmodified",...c.marks.map(m=>m.id)],pairs=states.flatMap((a,i)=>states.slice(i).map(b=>[a,b])),pair=pairs.find(([a,b])=>!c.chromatin.state_interactions.some(r=>(r.a===a&&r.b===b)||(r.a===b&&r.b===a)));if(!pair){toast("Every chromatin state pair already has a rule.");return;}c.chromatin.state_interactions.push({a:pair[0],b:pair[1],strength:.5,radius:1.5,min_separation:3});}
    else if(action==="remove-state-interaction")c.chromatin.state_interactions.splice(index,1);
    else if(action==="add-speckle"||action==="add-nucleolus") {const kind=action==="add-speckle"?"speckle":"nucleolus";c.bodies.push({id:uniqueId(kind),name:kind==="speckle"?`Speckle ${c.bodies.filter(b=>b.kind==="speckle").length+1}`:`Nucleolus ${c.bodies.filter(b=>b.kind==="nucleolus").length+1}`,kind,color:kind==="speckle"?"#b8a0cf":"#c8b385",center:[kind==="speckle"?3:-3,0,0],radius:kind==="speckle"?2:2.5,strength:1,chromatin_affinity:kind==="speckle"?0.3:-0.5,protein_affinity:Object.fromEntries(c.proteins.map(p=>[p.id,0])),enabled:true});}
    else if(action==="remove-body"){const [body]=c.bodies.splice(index,1);c.genes.forEach(g=>{if(g.body_responses)delete g.body_responses[body.id];});}
    else if(action==="remove-interaction")c.interactions.splice(index,1);
    else if(action==="add-nucleus"){const R=Number(c.simulation.confinement_radius)||20;c.nucleus={enabled:true,center:[-4.2*R,0,0],radii:[5*R,5*R,4*R],envelope_stiffness:30,tethers:[]};}
    else if(action==="remove-nucleus")delete c.nucleus;
    else if(action==="add-tether"){c.nucleus.tethers??=[];c.nucleus.tethers.push({id:uniqueId("tether"),name:`Tethering element ${c.nucleus.tethers.length+1}`,color:"#d6aa68",enabled:true,anchor:"envelope",start:0,end:Math.min(Math.max(1,Math.floor(c.chromatin.beads*.1)),c.chromatin.beads-1),attach_rate:.5,detach_rate:.01,capture_distance:1.5,spring_k:2,rest_length:.5,required_mark:null});}
    else if(action==="remove-tether")c.nucleus.tethers.splice(index,1);
    else if(action==="add-region") {if(!c.marks.length){toast("Add a mark type before defining its initial region.");return;}c.chromatin.initial_marks.push({start:0,end:Math.min(10,c.chromatin.beads-1),mark:c.marks[0].id});}
    else if(action==="remove-region")c.chromatin.initial_marks.splice(index,1);
    else if(action==="add-barrier")c.chromatin.barriers.push({bead:Math.floor(c.chromatin.beads/2),direction:"both"});
    else if(action==="remove-barrier")c.chromatin.barriers.splice(index,1);
    markDirty();renderEditor();
    const addedKind={"add-protein":"protein","add-mark":"mark","add-gene":"gene","add-speckle":"body","add-nucleolus":"body","add-enhancer":"regulatory","add-promoter":"regulatory","add-state-interaction":"state-interaction","add-tether":"tether"}[action];
    if(addedKind){const cards=$$("#editor-content details").filter(d=>d.dataset.card.startsWith(addedKind+"-")),last=cards.at(-1);if(last){last.open=true;requestAnimationFrame(()=>{last.scrollIntoView({block:"nearest",behavior:"smooth"});last.querySelector('input[type="text"]')?.focus({preventScroll:true});});}}
  }

  // Perspective canvas. Geometry comes only from the current server snapshot.
  function hexRGB(color) {let h=String(color||"#8ca28c").replace("#","");if(h.length===3)h=h.split("").map(x=>x+x).join("");const n=parseInt(h,16);return Number.isFinite(n)?[(n>>16)&255,(n>>8)&255,n&255]:[140,162,140];}
  function rgba(color,a){const rgb=hexRGB(color);return `rgba(${rgb[0]},${rgb[1]},${rgb[2]},${a})`;}
  function setupCanvas(canvas) {const rect=canvas.getBoundingClientRect(),dpr=Math.min(window.devicePixelRatio||1,2);const w=Math.max(1,rect.width),h=Math.max(1,rect.height);if(canvas.width!==Math.round(w*dpr)||canvas.height!==Math.round(h*dpr)){canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);}const ctx=canvas.getContext("2d");ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);return {ctx,w,h};}
  function makeProjector(w,h,radius) {
    const cy=Math.cos(S.camera.yaw),sy=Math.sin(S.camera.yaw),cp=Math.cos(S.camera.pitch),sp=Math.sin(S.camera.pitch);
    const base=Math.min(w*.43,h*.43)/Math.max(1,radius)*S.camera.zoom;
    return point=>{const [x,y,z]=point,xx=x*cy+z*sy,zz=-x*sy+z*cy,yy=y*cp-zz*sp,depth=y*sp+zz*cp,perspective=1/(1-depth/(radius*5));
      return {x:w*.5+xx*base*perspective,y:h*.49-yy*base*perspective,z:depth,scale:base*perspective};};
  }
  function drawScene() {
    if($("#scene-container").hidden)return;
    const {ctx,w,h}=setupCanvas($("#scene"));if(!S.snapshot||!S.applied)return;
    const c=S.applied, snap=S.snapshot,R=Number(c.simulation.confinement_radius)||12,env=snap.nucleus?.enabled?snap.nucleus:null,viewR=env&&S.view.nucleus&&S.view.fitNucleus?Math.max(R,...env.radii.map((r,k)=>r+Math.abs(env.center[k]))):R,proj=makeProjector(w,h,viewR);S.hitTargets=[];
    const nucleus=proj([0,0,0]),rr=Math.min(w*.43,h*.43)*S.camera.zoom*R/viewR;
    ctx.save();ctx.beginPath();ctx.ellipse(nucleus.x,nucleus.y+rr*.75,rr*.84,rr*.2,0,0,Math.PI*2);const shadow=ctx.createRadialGradient(nucleus.x,nucleus.y+rr*.75,0,nucleus.x,nucleus.y+rr*.75,rr);shadow.addColorStop(0,"rgba(83,127,169,.10)");shadow.addColorStop(1,"rgba(83,127,169,0)");ctx.fillStyle=shadow;ctx.fill();ctx.restore();
    const nucleusFill=ctx.createRadialGradient(nucleus.x-rr*.2,nucleus.y-rr*.2,0,nucleus.x,nucleus.y,rr);
    nucleusFill.addColorStop(0,"rgba(27,42,61,.58)");nucleusFill.addColorStop(1,"rgba(12,18,27,.45)");
    ctx.beginPath();ctx.arc(nucleus.x,nucleus.y,rr,0,Math.PI*2);ctx.fillStyle=nucleusFill;ctx.fill();ctx.strokeStyle="rgba(132,164,196,.22)";ctx.lineWidth=1;ctx.stroke();
    // The wire sphere makes the nuclear confinement explicit without obscuring particles.
    ctx.lineWidth=.7;ctx.strokeStyle="rgba(133,164,195,.14)";
    for(let latitude=-60;latitude<=60;latitude+=30){const a=latitude*Math.PI/180;ctx.beginPath();for(let j=0;j<=120;j++){const t=j*Math.PI*2/120,p=proj([R*Math.cos(a)*Math.cos(t),R*Math.sin(a),R*Math.cos(a)*Math.sin(t)]);if(j===0)ctx.moveTo(p.x,p.y);else ctx.lineTo(p.x,p.y);}ctx.stroke();}
    for(let meridian=0;meridian<180;meridian+=45){const a=meridian*Math.PI/180;ctx.beginPath();for(let j=0;j<=120;j++){const t=j*Math.PI*2/120,p=proj([R*Math.cos(t)*Math.cos(a),R*Math.sin(t),R*Math.cos(t)*Math.sin(a)]);if(j===0)ctx.moveTo(p.x,p.y);else ctx.lineTo(p.x,p.y);}ctx.stroke();}
    // Whole-nucleus envelope: a separate soft ellipsoid around the modeled territory.
    if(env&&S.view.nucleus){const [ax,ay,az]=env.radii,[cx,cy0,cz]=env.center,ring=(fn)=>{ctx.beginPath();for(let j=0;j<=160;j++){const t=j*Math.PI*2/160,p=proj(fn(t));if(j===0)ctx.moveTo(p.x,p.y);else ctx.lineTo(p.x,p.y);}ctx.stroke();};ctx.strokeStyle="rgba(214,170,104,.42)";ctx.lineWidth=1.4;ring(t=>[cx+ax*Math.cos(t),cy0+ay*Math.sin(t),cz]);ring(t=>[cx+ax*Math.cos(t),cy0,cz+az*Math.sin(t)]);ring(t=>[cx,cy0+ay*Math.cos(t),cz+az*Math.sin(t)]);ctx.strokeStyle="rgba(214,170,104,.16)";ctx.lineWidth=.7;for(const lat of [-.5,.5]){const s=Math.sqrt(1-lat*lat);ring(t=>[cx+ax*s*Math.cos(t),cy0+ay*lat,cz+az*s*Math.sin(t)]);}if(S.view.labels){const p=proj([cx,cy0+ay,cz]);ctx.font="500 9px -apple-system,sans-serif";ctx.fillStyle="rgba(226,186,120,.9)";ctx.textAlign="center";ctx.fillText("NUCLEAR ENVELOPE",p.x,p.y-8);}}
    const vacant=new Set();if(S.view.occupancy&&snap.nucleosomes?.enabled)for(const [a,b,v] of snap.genome_tracks?.nucleosome_occupancy_runs||[])if(!v)for(let i=a;i<b;i++)vacant.add(i);
    const proteinById=new Map((snap.proteins||[]).map(p=>[p.id,p])),speciesById=new Map(c.proteins.map(p=>[p.id,p]));
    const positions=snap.chromatin||[],projected=positions.map(proj),detailScale=positions.length>700?Math.max(.45,Math.sqrt(600/positions.length)):1,beadColors=new Array(positions.length).fill("#738c9e");
    for(const mark of c.marks){const flags=snap.marks?.[mark.id]||[];for(let i=0;i<positions.length;i++)if(flags[i])beadColors[i]=mark.color;}
    const selectedGene=S.selected?.type==="gene"?c.genes.find(g=>g.id===S.selected.id):S.selected?.type==="regulatory"?c.regulatory_regions.find(r=>r.id===S.selected.id):null;
    const geneHighlighted=i=>selectedGene&&i>=selectedGene.start&&i<=selectedGene.end;
    const objects=[];
    for(let i=0;i<positions.length-1;i++){const a=projected[i],b=projected[i+1];objects.push({z:(a.z+b.z)/2,draw:()=>{ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.strokeStyle=geneHighlighted(i)?"#f08090":rgba(beadColors[i],.67);ctx.lineWidth=geneHighlighted(i)?2.5:Math.max(.65,Math.min(3,a.scale*.14)*detailScale);ctx.stroke();}});}
    if(S.view.proteins)for(const chain of snap.protein_chains||[]){for(let j=0;j<chain.length-1;j++){const aa=proteinById.get(chain[j]),bb=proteinById.get(chain[j+1]);if(!aa||!bb)continue;const a=proj(aa.position),b=proj(bb.position),color=speciesById.get(aa.species)?.color||"#b4a3c7";objects.push({z:(a.z+b.z)/2,draw:()=>{ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.strokeStyle=rgba(color,.55);ctx.lineWidth=1.1;ctx.stroke();}});}}
    if(S.view.bonds)for(const bond of snap.bonds||[]){const a=proj(bond.a),b=proj(bond.b),isCoh=String(bond.kind).includes("cohes"),isTether=bond.kind==="tether";if(isTether&&!S.view.nucleus)continue;objects.push({z:(a.z+b.z)/2,draw:()=>{ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.strokeStyle=isTether?"rgba(226,186,120,.85)":isCoh?"rgba(195,146,68,.7)":"rgba(72,147,126,.42)";ctx.lineWidth=isTether?1.4:isCoh?1.8:.8;ctx.setLineDash(isCoh||isTether?[]:[3,2]);ctx.stroke();ctx.setLineDash([]);if(isTether){ctx.fillStyle="rgba(226,186,120,.95)";ctx.fillRect(b.x-2.5,b.y-2.5,5,5);}}});}
    c.bodies.filter(b=>b.enabled).forEach((body,index)=>{const p=proj(body.center),r=Number(body.radius)*p.scale;objects.push({z:p.z-r*.02,draw:()=>{const grad=ctx.createRadialGradient(p.x-r*.26,p.y-r*.26,0,p.x,p.y,r);grad.addColorStop(0,rgba(body.color,.20));grad.addColorStop(.65,rgba(body.color,.13));grad.addColorStop(1,rgba(body.color,.07));ctx.beginPath();ctx.arc(p.x,p.y,r,0,Math.PI*2);ctx.fillStyle=grad;ctx.fill();ctx.strokeStyle=rgba(body.color,.34);ctx.lineWidth=1;ctx.stroke();ctx.beginPath();ctx.ellipse(p.x,p.y,r*.84,r*.32,-.2,0,Math.PI*2);ctx.strokeStyle=rgba(body.color,.12);ctx.stroke();if(S.view.labels){ctx.fillStyle=rgba(body.color,.94);ctx.font="500 9px -apple-system, sans-serif";ctx.textAlign="center";ctx.fillText(body.name,p.x,p.y+r+13);ctx.font="7px -apple-system, sans-serif";ctx.fillStyle="#a1a994";ctx.fillText("IMPOSED FIELD",p.x,p.y+r+24);}S.hitTargets.push({x:p.x,y:p.y,r:Math.max(8,r*.7),z:p.z,type:"body",id:body.id});}});});
    positions.forEach((position,i)=>{const p=projected[i],r=Math.max(1,Math.max(2,Math.min(5.6,p.scale*.17))*detailScale),selected=(S.selected?.type==="chromatin"&&S.selected.index===i)||(S.selected?.type==="loop"&&(S.selected.left===i||S.selected.right===i))||geneHighlighted(i);objects.push({z:p.z,draw:()=>{if(selected){ctx.beginPath();ctx.arc(p.x,p.y,r+5,0,Math.PI*2);ctx.strokeStyle="#9fcaec";ctx.lineWidth=1.5;ctx.stroke();}ctx.beginPath();ctx.arc(p.x,p.y,r,0,Math.PI*2);if(vacant.has(i)){ctx.fillStyle="rgba(10,16,24,.85)";ctx.fill();ctx.strokeStyle="rgba(221,213,230,.9)";ctx.lineWidth=1;ctx.stroke();}else{ctx.fillStyle=rgba(beadColors[i],.75+(p.z+R)/(R*2)*.2);ctx.fill();if(positions.length<700||selected){ctx.strokeStyle="rgba(255,255,255,.65)";ctx.lineWidth=.6;ctx.stroke();}}S.hitTargets.push({x:p.x,y:p.y,r:Math.max(7,r+2),z:p.z,type:"chromatin",index:i});}});});
    if(S.view.proteins)(snap.proteins||[]).forEach((particle,i)=>{const species=speciesById.get(particle.species);if(!species)return;const p=proj(particle.position),r=Math.max(1.8,Math.min(6,p.scale*(Number(species.radius)||.25)*.48)),selected=S.selected?.type==="protein"&&S.selected.index===i;objects.push({z:p.z,draw:()=>{if(selected){ctx.beginPath();ctx.arc(p.x,p.y,r+5,0,Math.PI*2);ctx.strokeStyle=species.color;ctx.lineWidth=1.5;ctx.stroke();}ctx.beginPath();ctx.arc(p.x,p.y,r,0,Math.PI*2);ctx.fillStyle=rgba(species.color,particle.sticker ? .88 : .5);ctx.fill();if(particle.sticker){ctx.strokeStyle=rgba(species.color,.9);ctx.lineWidth=.6;ctx.stroke();}S.hitTargets.push({x:p.x,y:p.y,r:Math.max(6,r+2),z:p.z,type:"protein",index:i});}});});
    if(S.view.rna&&snap.transcription?.enabled)for(const transcript of snap.transcription.transcripts||[]){if(!Array.isArray(transcript.position))continue;const p=proj(transcript.position),nascent=transcript.state==="nascent",color=window.ChromatinTranscriptHelpers.transcriptColor(transcript),selected=S.selected?.type==="transcript"&&S.selected.id===transcript.id,r=nascent?4:3;objects.push({z:p.z,draw:()=>{if(selected){ctx.beginPath();ctx.arc(p.x,p.y,r+5,0,Math.PI*2);ctx.strokeStyle=color;ctx.lineWidth=1.5;ctx.stroke();}ctx.fillStyle=color;ctx.strokeStyle="#fff";ctx.lineWidth=.7;ctx.beginPath();if(nascent){ctx.moveTo(p.x,p.y-r-1);ctx.lineTo(p.x+r,p.y);ctx.lineTo(p.x,p.y+r+1);ctx.lineTo(p.x-r,p.y);ctx.closePath();}else ctx.arc(p.x,p.y,r,0,Math.PI*2);ctx.fill();ctx.stroke();S.hitTargets.push({x:p.x,y:p.y,r:Math.max(6,r+2),z:p.z,type:"transcript",id:transcript.id});}});}
    objects.sort((a,b)=>a.z-b.z).forEach(o=>o.draw());
    // Locus labels and current extruder endpoints remain readable above the particle layer.
    if(S.view.labels)c.genes.forEach(g=>{const index=Math.round((g.start+g.end)/2);if(!positions[index])return;const p=projected[index],state=snap.genes?.find(s=>s.id===g.id),label=g.name;ctx.font="500 9px -apple-system,sans-serif";const tw=ctx.measureText(label).width;ctx.fillStyle=state?.on?"rgba(66,29,42,.96)":"rgba(19,28,40,.96)";roundRect(ctx,p.x+10,p.y-19,tw+15,18,4);ctx.fill();ctx.strokeStyle=state?.on?"rgba(240,128,144,.5)":"rgba(133,164,195,.4)";ctx.lineWidth=.7;ctx.stroke();ctx.fillStyle=state?.on?"#ffabb6":"#c1cedd";ctx.textAlign="left";ctx.fillText(label,p.x+17,p.y-7);ctx.beginPath();ctx.moveTo(p.x+2,p.y-2);ctx.lineTo(p.x+10,p.y-10);ctx.strokeStyle="#748ba3";ctx.stroke();});
    const cohesinRing=(snap.cohesin||[]).length>40?2.6:5;if(S.view.bonds)(snap.cohesin||[]).forEach(([a,b])=>{if(!positions[a]||!positions[b])return;for(const idx of [a,b]){const p=projected[idx];ctx.beginPath();ctx.arc(p.x,p.y,cohesinRing,0,Math.PI*2);ctx.strokeStyle="#c29452";ctx.lineWidth=1.5;ctx.stroke();}});
    // A scale reference is geometry-derived, not a biological calibration.
    const unit=Math.pow(10,Math.floor(Math.log10(R/3))),len=unit*nucleus.scale;
    if(len<100&&len>5){ctx.strokeStyle="#a7b59c";ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(w-80,h-55);ctx.lineTo(w-80+len,h-55);ctx.moveTo(w-80,h-58);ctx.lineTo(w-80,h-52);ctx.moveTo(w-80+len,h-58);ctx.lineTo(w-80+len,h-52);ctx.stroke();ctx.font="9px -apple-system,sans-serif";ctx.textAlign="left";ctx.fillStyle="#a1ae95";ctx.fillText(`${pretty(unit)} length units`,w-80,h-63);}
  }
  function roundRect(ctx,x,y,w,h,r){ctx.beginPath();ctx.moveTo(x+r,y);ctx.arcTo(x+w,y,x+w,y+h,r);ctx.arcTo(x+w,y+h,x,y+h,r);ctx.arcTo(x,y+h,x,y,r);ctx.arcTo(x,y,x+w,y,r);ctx.closePath();}
  function queueRender(){if(S.renderPending)return;S.renderPending=true;requestAnimationFrame(()=>{S.renderPending=false;drawScene();if(S.chartsDirty){drawCharts();S.chartsDirty=false;}});}
  function drawCharts() {
    chart($("#chart-fractions"),[{key:"bound_fraction",color:"#79c9bd"},{key:"marked_fraction",color:"#b19cda"}],{fixedMax:1});
    chart($("#chart-structure"),[{key:"largest_cluster",color:"#c78b44"},{key:"rg",color:"#88b4dc",right:true}],{});
    chart($("#chart-rna"),[{key:"nascent",color:"#a890c0"},{key:"unprocessed",color:"#d3a45e"},{key:"mature",color:"#79c9bd"}],{});
  }
  function chart(canvas,series,options) {
    const {ctx,w,h}=setupCanvas(canvas),left=25,right=series.some(s=>s.right)?25:8,top=13,bottom=20,pw=w-left-right,ph=h-top-bottom;
    ctx.font="9px -apple-system,sans-serif";ctx.textAlign="right";ctx.fillStyle="#a1b2c5";ctx.lineWidth=.6;
    const history=S.history;
    const maxFor=(rightAxis)=>Math.max(options.fixedMax||1,...history.flatMap(row=>series.filter(s=>!!s.right===rightAxis).map(s=>Number(row[s.key])||0)))*(!options.fixedMax?1.1:1);
    const max=maxFor(false),rightMax=maxFor(true),end=history.at(-1)?.time||1,start=history[0]?.time||0;
    for(let i=0;i<3;i++){const y=top+i*ph/2;ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(w-right,y);ctx.strokeStyle="#2c3644";ctx.stroke();ctx.fillStyle="#a1b2c5";ctx.textAlign="right";ctx.fillText(pretty(max*(1-i/2),max<=2?1:0),left-5,y+2);if(series.some(s=>s.right)){ctx.textAlign="left";ctx.fillStyle="#9ab8d3";ctx.fillText(pretty(rightMax*(1-i/2),1),w-right+4,y+2);}}
    ctx.textAlign="left";ctx.fillStyle="#a1b2c5";const timeDigits=end-start<.01?5:end-start<.1?3:end-start<1?2:1;ctx.fillText(pretty(start,timeDigits),left,h-5);ctx.textAlign="right";ctx.fillText(`${pretty(end,timeDigits)} t`,w-right,h-5);
    if(history.length<2){ctx.textAlign="center";ctx.fillStyle="#a1b2c5";ctx.font="9px -apple-system,sans-serif";ctx.fillText("Run to collect a trajectory",left+pw/2,top+ph*.55);return;}
    for(const s of series){ctx.beginPath();history.forEach((row,i)=>{const x=left+(row.time-start)/Math.max(1e-9,end-start)*pw,y=top+ph-(Number(row[s.key])||0)/(s.right?rightMax:max)*ph;if(i===0)ctx.moveTo(x,y);else ctx.lineTo(x,y);});ctx.strokeStyle=s.color;ctx.lineWidth=1.45;ctx.lineJoin="round";ctx.stroke();}
  }
  function renderLegend(){if(!S.applied)return;const c=S.applied;const items=[["#738c9e","Chromatin",false],...c.marks.map(m=>[m.color,m.name,false]),...c.proteins.map(p=>[p.color,p.name,false]),...c.bodies.filter(b=>b.enabled).map(b=>[b.color,`${b.name} · field`,true]),...(c.transcription?.enabled?[["#8c70b4","Polymerase / nascent RNA",false],["#b98549","Paused polymerase",false],["#839399","Stalled polymerase",false],["#d09d4d","Released RNA",false],["#138d91","Processed RNA",false]]:[])];$("#legend").innerHTML=items.map(([color,name,landmark])=>`<span class="legend-item"><i class="${landmark?"landmark":""}" style="background:${esc(color)};--legend-color:${esc(color)}"></i>${esc(name)}</span>`).join("");}
  function prop(label,value){return `<div class="inspector-property"><span>${esc(label)}</span><strong>${esc(value)}</strong></div>`;}
  function coords(p){return `<div class="inspector-coords">x ${pretty(p[0],2)} &nbsp; y ${pretty(p[1],2)} &nbsp; z ${pretty(p[2],2)}</div>`;}
  function renderInspector(){
    const s=S.selected,c=S.applied,snap=S.snapshot;
    $("#clear-selection").hidden=!s;$("#selection-chip").hidden=!s;
    if(!s||!c||!snap){$("#inspector-content").innerHTML='<div class="inspector-icon">⌖</div><h3>Explore the scene</h3><p>Select a chromatin bead, protein, or landmark to inspect its state.</p><div class="inspector-note">Motion and binding are explicit. Colored nuclear bodies are imposed affinity fields.</div>';return;}
    let html="",label="";
    if(s.type==="chromatin"){
      const p=snap.chromatin[s.index];if(!p){S.selected=null;renderInspector();return;}
      label=`Chromatin bead ${s.index}`;const marks=c.marks.filter(m=>snap.marks?.[m.id]?.[s.index]),genes=c.genes.filter(g=>s.index>=g.start&&s.index<=g.end);
      html=`<div class="inspector-icon">⌁</div><h3>${label}</h3>`+prop("Genomic start",`${pretty(s.index*c.chromatin.bp_per_bead,0)} bp`)+prop("Inclusive end",`${pretty((s.index+1)*c.chromatin.bp_per_bead-1,0)} bp`)+prop("Bead resolution",genomicSize(c.chromatin.bp_per_bead))+prop("Marks",marks.length?marks.map(m=>m.name).join(", "):"Unmodified")+prop("Gene coverage",genes.map(g=>g.name).join(", ")||"None")+coords(p);
      const distances=c.bodies.filter(b=>b.enabled).map(b=>({name:b.name,d:Math.sqrt(p.reduce((sum,v,i)=>sum+(v-b.center[i])**2,0))-b.radius})).sort((a,b)=>a.d-b.d);
      if(distances.length)html+=prop("Nearest field surface",`${pretty(distances[0].d,2)} · ${distances[0].name}`);
      html+='<p class="inspector-note">Genomic coordinates are zero-based model coordinates. Negative surface distance indicates a bead inside an imposed field.</p>';
    }else if(s.type==="protein"){
      const p=snap.proteins[s.index],species=p&&c.proteins.find(x=>x.id===p.species);if(!p||!species){S.selected=null;renderInspector();return;}
      label=species.name;html=`<div class="inspector-icon" style="color:${esc(species.color)}">●</div><h3>${esc(species.name)}</h3>`+prop("Configured role",proteinRoles(c,species.id).join(" · ")||"Other protein species")+prop("Molecule",p.molecule)+prop("Particle",p.id??s.index)+prop("Representation",species.representation==="idr"?"Coarse IDR chain":"Single bead")+prop("Binding site",p.sticker?"Sticker":"Non-sticker")+prop("Free diffusion D",pretty(species.diffusion,3))+prop("Molecular valence",species.valence)+coords(p.position)+'<p class="inspector-note">Observed motion reflects links, interactions, and confinement as well as free diffusion.</p>';
      const molecule=snap.interaction_state?.molecules?.find(m=>m.id===p.molecule);
      if(molecule){const partners=(snap.interaction_state.molecule_links||[]).filter(link=>link.a===molecule.id||link.b===molecule.id);html+=prop("DNA attachment",molecule.binding_status==="direct"?"Directly bound":molecule.binding_status==="indirect"?"Indirectly tethered":"Not tethered")+prop("Linked protein molecules",partners.length)+prop("Valence used",`${molecule.used_valence} / ${molecule.valence}`);}
      html+=`<button class="secondary-button" data-inspect-interactions="${esc(species.id)}">View protein interactions</button>`;
    }else if(s.type==="transcript"){
      const transcript=snap.transcription?.transcripts?.find(t=>t.id===s.id),gene=transcript&&c.genes.find(g=>g.id===transcript.gene_id);if(!transcript){S.selected=null;renderInspector();return;}
      label=`${transcript.id} · ${transcript.state}`;html=`<div class="inspector-icon" style="color:#8c70b4">≋</div><h3>${esc(transcript.id)}</h3>`+prop("Gene",gene?.name||transcript.gene_id)+prop("RNA state",transcript.state)+prop("Length",`${pretty(transcript.length_nt,0)} nt`)+prop("Polymerase",transcript.polymerase_id)+prop("Polymerase state",transcript.state==="nascent"?window.ChromatinTranscriptHelpers.polymeraseMotion(transcript):"Released (origin ID)")+(transcript.state==="nascent"?prop("Genomic position",`${pretty(transcript.bead_position*c.chromatin.bp_per_bead,0)} bp`)+prop("Current elongation",`${pretty(transcript.elongation_rate*c.chromatin.bp_per_bead,2)} nt / time`):"")+coords(transcript.position)+`<button class="secondary-button" data-inspect-transcript="${esc(transcript.id)}">View elongation & RNA</button>`+'<p class="inspector-note">RNA is an individual passive point. Nascent RNA follows its polymerase; released RNA diffuses. No RNA chain geometry is implied.</p>';
    }else if(s.type==="body"){
      const b=c.bodies.find(x=>x.id===s.id);if(!b){S.selected=null;renderInspector();return;}
      label=b.name;html=`<div class="inspector-icon" style="color:${esc(b.color)}">◉</div><h3>${esc(b.name)}</h3>`+prop("Model",`Imposed ${b.kind} field`)+prop("Radius",pretty(b.radius,2))+prop("Field strength",pretty(b.strength,2))+prop("Chromatin affinity",pretty(b.chromatin_affinity,2))+coords(b.center)+'<p class="inspector-note">This is a fixed spatial field. Its geometry does not emerge from protein assembly.</p>';
    }else if(s.type==="loop"){
      const active=(snap.cohesin||[]).some(([a,b])=>a===s.left&&b===s.right);
      label=`Anchors ${s.left} ↔ ${s.right}`;html=`<div class="inspector-icon">⋈</div><h3>${esc(label)}</h3>`+prop("Anchor pair",active?"Current loop":"No longer active")+prop("Left genomic start",`${pretty(s.left*c.chromatin.bp_per_bead,0)} bp`)+prop("Right inclusive end",`${pretty((s.right+1)*c.chromatin.bp_per_bead-1,0)} bp`)+prop("Covered span",genomicSize((s.right-s.left+1)*c.chromatin.bp_per_bead))+'<p class="inspector-note">Anchor positions identify this pair, not a persistent motor lineage. Both anchors are highlighted in the 3D view.</p>';
    }else if(s.type==="regulatory"){
      const region=c.regulatory_regions.find(r=>r.id===s.id),live=snap.regulation?.regions?.find(r=>r.id===s.id);if(!region){S.selected=null;renderInspector();return;}
      label=region.name;html=`<div class="inspector-icon" style="color:${esc(region.color)}">⇢</div><h3>${esc(region.name)}</h3>`+prop("Element",region.kind)+prop("Target gene",c.genes.find(g=>g.id===region.target_gene)?.name||region.target_gene)+prop("Genomic interval",`${pretty(region.start*c.chromatin.bp_per_bead,0)}–${pretty((region.end+1)*c.chromatin.bp_per_bead,0)} bp`)+prop("Bound TF molecules",live?live.tf_count:"Unavailable")+prop("Recruited coactivators",live?live.recruited_count:"Unavailable")+prop("TF engaged",live?live.active?"Yes":"No":"Unavailable")+prop("Contact gate",live?live.contact?"In contact":"Out of contact":"Unavailable")+`<button class="secondary-button" data-inspect-regulation="${esc(region.id)}">View regulatory mechanism</button>`;
    }else if(s.type==="gene"){
      const g=c.genes.find(x=>x.id===s.id),state=snap.genes.find(x=>x.id===s.id);if(!g||!state){S.selected=null;renderInspector();return;}
      label=g.name;html=`<div class="inspector-icon">≋</div><h3>${esc(g.name)}</h3>`+prop("Promoter",`${state.on?"ON":"OFF"} · stochastic switching`)+prop("RNA model",c.transcription?.enabled?"Individual polymerases and RNA":"Aggregate RNA counts")+prop("Bead interval",`${g.start}–${g.end}`)+prop("Genomic start",`${pretty(g.start*c.chromatin.bp_per_bead,0)} bp`)+prop("Inclusive end",`${pretty((g.end+1)*c.chromatin.bp_per_bead-1,0)} bp`)+prop("Gene span",genomicSize((g.end-g.start+1)*c.chromatin.bp_per_bead))+prop("Elongation",`${pretty(g.elongation_rate*c.chromatin.bp_per_bead/1000,2)} kb / reduced time`)+prop("Nascent polymerases",state.nascent)+prop("Unprocessed RNA",state.unprocessed)+prop("Processed RNA",state.mature)+prop("Total initiated",state.initiated)+prop("Total completed",state.completed)+prop("Speckle contact",state.speckle_contact?"Yes":"No")+prop("Active stiffness gain",`${pretty(g.persistence_multiplier,2)}×`)+'<p class="inspector-note">Speckle processing is a downstream RNA effect. Separate optional body-response hypotheses can change promoter activation and initiation; Regulation shows their current factors.</p>';
    }
    if(["chromatin","gene","loop","regulatory","transcript"].includes(s.type))html+=genomeView?.inspectorHTML(s)||"";
    $("#inspector-content").innerHTML=html;$("#selection-chip").textContent=label;
  }
  function renderGenes(){const genes=S.snapshot?.genes||[];$("#gene-count").textContent=`${genes.length} ${genes.length===1?"locus":"loci"}`;$("#gene-activity").innerHTML=genes.length?genes.map(g=>`<div class="gene-row"><div class="gene-name-row"><button data-select-gene="${esc(g.id)}" title="Inspect ${esc(g.name)}">${esc(g.name)}</button><span class="gene-state ${g.on?"on":""}">${g.on?"ON":"OFF"}</span></div><div class="gene-numbers"><span><strong>${pretty(g.nascent,0)}</strong>nascent</span><span><strong>${pretty(g.unprocessed,0)}</strong>unproc.</span><span><strong>${pretty(g.mature,0)}</strong>processed</span></div></div>`).join(""):'<p class="small-muted">No genes configured. Add a locus in the model builder.</p>';}
  function renderEvents(){const events=S.events.slice(-60).reverse();$("#event-count").textContent=S.events.length;$("#events").innerHTML=events.length?events.map(e=>{const type=String(e.type||"event"),group=/bind|link|bond/.test(type)?"binding":/trans|init|rna|promot|elong/.test(type)?"transcription":/cohes|extrud/.test(type)?"cohesin":"mark",symbol=group==="binding"?"⋈":group==="transcription"?"≋":group==="cohesin"?"↔":"✦";return `<div class="event-item"><span class="event-symbol ${group}">${symbol}</span><div class="event-content"><div class="event-topline"><strong>${esc(type.replace(/_/g," "))}</strong><time>${pretty(e.time,2)} t</time></div><p>${esc(typeof e.detail==="string"?e.detail:JSON.stringify(e.detail))}</p></div></div>`;}).join(""):'<p class="empty-events">Binding, catalysis, and transcription events will appear as the simulation runs.</p>';}

  function selectScene(selection,reveal=false){S.selected=selection;genomeView?.setSelection(selection,{reveal});regulationView?.setSelection(selection);renderInspector();queueRender();}
  function showSceneView(view){S.sceneView=view;const descriptions={"3d":"Follow the actual chromatin and proteins in 3D. Select a molecule or locus to inspect it.","1d":"Read the applied genome from left to right: genes, loop anchors, marks, access and binding preferences.",interactions:"Inspect which proteins can bind, which are linked now, and whether an assembly touches chromatin.",regulation:"Connect local chromatin and regulatory inputs to promoter switching, elongation and RNA output."};$("#view-explanation").textContent=descriptions[view]||descriptions["3d"];for(const [name,container] of [["3d","scene-container"],["1d","genome-container"],["interactions","interactions-container"],["regulation","regulation-container"]]){$(`#${container}`).hidden=view!==name;$(`#view-${name}`).classList.toggle("active",view===name);$(`#view-${name}`).setAttribute("aria-pressed",String(view===name));}if(view==="1d")genomeView?.schedule();else queueRender();}
  genomeView=new window.ChromatinGenomeView($("#genome-container"),{onSelect:selection=>selectScene(selection),onRun:safe(run),onStep:safe(stepOnce)});
  interactionView=new window.ChromatinInteractionView($("#interactions-container"),{estimate:request=>api("/api/kinetics/estimate",request)});
  regulationView=new window.ChromatinRegulationView($("#regulation-container"),{onSelect:selection=>selectScene(selection),onNavigate:(view,selection)=>{selectScene(selection,true);if(view==="1d")genomeView?.action("selected-bin");showSceneView(view);$("#mobile-workspace-nav [data-mobile-tab='scene']").click();},onConfigure:()=>{$("#editor-nav [data-tab='regulation']").click();$("#mobile-workspace-nav [data-mobile-tab='editor']").click();}});
  $("#view-3d").addEventListener("click",()=>showSceneView("3d"));$("#view-1d").addEventListener("click",()=>showSceneView("1d"));$("#view-interactions").addEventListener("click",()=>showSceneView("interactions"));$("#view-regulation").addEventListener("click",()=>showSceneView("regulation"));
  function download(data,filename){const blob=new Blob([JSON.stringify(data,null,2)],{type:"application/json"}),url=URL.createObjectURL(blob),a=document.createElement("a");a.href=url;a.download=filename;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  function fileName(suffix){return `${(S.config?.name||"chromatin-studio").toLowerCase().replace(/[^a-z0-9]+/g,"-").replace(/^-|-$/g,"")||"chromatin-studio"}-${suffix}.json`;}
  $("#mobile-workspace-nav").addEventListener("click",e=>{const b=e.target.closest("[data-mobile-tab]");if(!b)return;if(b.dataset.mobileTab!=="editor"&&S.builderExpanded)setBuilderExpanded(false);if(["learn","calibration","studies"].includes(S.tab)&&b.dataset.mobileTab!=="editor"){S.tab="system";$$("#editor-nav button").forEach(x=>x.classList.toggle("active",x.dataset.tab===S.tab));renderEditor();}document.body.dataset.mobileView=b.dataset.mobileTab;$$("#mobile-workspace-nav button").forEach(x=>{x.classList.toggle("active",x===b);x.setAttribute("aria-pressed",String(x===b));});queueRender();});
  function focusHeading(selector){const heading=$(selector);if(heading){heading.setAttribute("tabindex","-1");heading.focus({preventScroll:true});}}
  function updateBuilderLayout(){
    const expanded=Boolean(S.builderExpanded),button=$("#expand-builder-btn");
    document.body.classList.toggle("builder-workspace",expanded);
    button.hidden=["cluster","learn","calibration","studies"].includes(S.tab);
    button.textContent=expanded?"← Return to simulation":"Expand builder ↗";
    button.setAttribute("aria-expanded",String(expanded));
  }
  function setBuilderExpanded(expanded){S.builderExpanded=expanded;updateBuilderLayout();S.chartsDirty=true;queueRender();window.scrollTo({top:0,behavior:"instant"});}
  $("#expand-builder-btn").addEventListener("click",()=>setBuilderExpanded(!S.builderExpanded));
  function navigateEditor(section){if(section==="simulation"){setBuilderExpanded(false);navigateEditor("system");$("#mobile-workspace-nav [data-mobile-tab='scene']").click();focusHeading("#experiment-title");return;}const button=$(`#editor-nav [data-tab="${section}"]`);if(!button)return;button.click();$("#mobile-workspace-nav [data-mobile-tab='editor']").click();focusHeading("#editor-content h2");}
  function openLearning(topic="overview"){S.learnTopic=topic;navigateEditor("learn");window.scrollTo({top:0,behavior:"instant"});}
  $("#learn-btn").addEventListener("click",()=>openLearning());
  $("#calibration-status").addEventListener("click",()=>navigateEditor("calibration"));
  $("#editor-content").addEventListener("click",e=>{const link=e.target.closest("[data-open-lesson]");if(link)openLearning(link.dataset.openLesson);});
  $(".intro-steps").addEventListener("click",e=>{const button=e.target.closest("[data-workflow]");if(!button)return;if(button.dataset.workflow==="build"){navigateEditor("system");if(window.innerWidth>700)setBuilderExpanded(true);}if(button.dataset.workflow==="test")navigateEditor("studies");if(button.dataset.workflow==="observe"){navigateEditor("system");$("#mobile-workspace-nav [data-mobile-tab='scene']").click();$("#view-3d").click();$(".simulation-panel").scrollIntoView({block:"start",behavior:"instant"});}});
  $("#view-explain-btn").addEventListener("click",()=>openLearning(S.sceneView==="regulation"?"transcription":S.sceneView==="interactions"?"binding":S.sceneView==="1d"?"chromatin":"motion"));
  $("#editor-nav").addEventListener("click",e=>{const b=e.target.closest("[data-tab]");if(!b)return;S.tab=b.dataset.tab;$$("#editor-nav button").forEach(x=>x.classList.toggle("active",x===b));$("#editor-content").innerHTML="";$("#editor-content").scrollTop=0;renderEditor();window.scrollTo({top:0,behavior:"instant"});});
  $("#editor-content").addEventListener("input",e=>{
    const el=e.target;
    if(!el.matches('input[type="text"],input[type="number"]'))return;
    // Store the draft immediately without replacing the focused input. A partially
    // typed numeric value remains invalid (null), never the stale applied value.
    const numeric=el.type==="number";
    const value=numeric?(el.value.trim()!==""&&Number.isFinite(Number(el.value))?Number(el.value):null):el.value;
    if(el.dataset.path){setPath(JSON.parse(el.dataset.path),value);markDirty();}
    if(el.dataset.pairA){
      const a=el.dataset.pairA,b=el.dataset.pairB;
      let pair=S.config.interactions.find(x=>(x.a===a&&x.b===b)||(x.a===b&&x.b===a));
      if(pair)pair.strength=value;
      else if(value!==0)S.config.interactions.push({a,b,strength:value,on_rate:0.3,off_rate:0.1,radius:1.2});
      $$("#editor-content [data-pair-a]").forEach(other=>{if(other!==el&&((other.dataset.pairA===a&&other.dataset.pairB===b)||(other.dataset.pairA===b&&other.dataset.pairB===a)))other.value=value??"";});
      markDirty();
    }
  });
  $("#editor-content").addEventListener("change",e=>{
    const el=e.target;
    if(el.dataset.path){const path=JSON.parse(el.dataset.path);let value=el.type==="checkbox"?el.checked:el.type==="number"?(el.value.trim()===""?null:Number(el.value)):el.value;if(el.type==="number"&&!Number.isFinite(value)){setPath(path,null);markDirty();return;}setPath(path,value);markDirty();if(["representation","name","color"].includes(path.at(-1)))renderEditor();}
    if(el.dataset.pairA){const a=el.dataset.pairA,b=el.dataset.pairB,strength=el.value.trim()===""?null:Number(el.value);let pair=S.config.interactions.find(x=>(x.a===a&&x.b===b)||(x.a===b&&x.b===a));if(!Number.isFinite(strength)||strength<0){showError("Interaction strengths must be nonnegative numbers.");return;}if(pair)pair.strength=strength;else if(strength>0)S.config.interactions.push({a,b,strength,on_rate:0.3,off_rate:0.1,radius:1.2});markDirty();renderEditor();}
  });
  $("#editor-content").addEventListener("click",e=>{const b=e.target.closest("[data-action]");if(!b)return;e.preventDefault();e.stopPropagation();handleAction(b.dataset.action,Number(b.dataset.index));});
  $("#run-btn").addEventListener("click",safe(run));
  $("#step-btn").addEventListener("click",safe(stepOnce));
  $("#apply-btn").addEventListener("click",safe(()=>reset(S.config)));
  $("#reset-btn").addEventListener("click",safe(()=>{S.backend=S.appliedBackend||S.backend;$("#engine-select").value=S.backend;return reset(S.applied);}));
  $("#dismiss-error").addEventListener("click",()=>$("#error-banner").hidden=true);
  $("#clear-selection").addEventListener("click",()=>selectScene(null));
  $("#inspector-content").addEventListener("click",e=>{const button=e.target.closest("[data-inspect-interactions]");if(button){interactionView?.focusSpecies(button.dataset.inspectInteractions);showSceneView("interactions");$("#mobile-workspace-nav [data-mobile-tab='scene']").click();}});
  $("#inspector-content").addEventListener("click",e=>{const button=e.target.closest("[data-inspect-regulation]");if(button){selectScene({type:"regulatory",id:button.dataset.inspectRegulation});showSceneView("regulation");$("#mobile-workspace-nav [data-mobile-tab='scene']").click();}});
  $("#inspector-content").addEventListener("click",e=>{const button=e.target.closest("[data-inspect-transcript]");if(button){selectScene({type:"transcript",id:button.dataset.inspectTranscript});showSceneView("regulation");$("#mobile-workspace-nav [data-mobile-tab='scene']").click();}});
  $("#gene-activity").addEventListener("click",e=>{const b=e.target.closest("[data-select-gene]");if(b){selectScene({type:"gene",id:b.dataset.selectGene},true);}});
  $("#preset-select").addEventListener("change",safe(async e=>{if(!e.target.value)return;pause();const id=e.target.value,config=await api(`/api/presets/${encodeURIComponent(id)}`);S.presetId=id;S.config=normalizeConfig(config);S.dirty=true;S.tab="system";$$("#editor-nav button").forEach(b=>b.classList.toggle("active",b.dataset.tab===S.tab));await reset(S.config,{resetCamera:true});toast("Preset loaded. Run the simulation to explore its dynamics.");}));
  $("#engine-select").addEventListener("change",e=>{S.backend=e.target.value;markDirty({preservePreset:true});toast("Apply & restart to initialize the selected engine.");});
  $("#import-btn").addEventListener("click",()=>$("#import-file").click());
  $("#import-file").addEventListener("change",safe(async e=>{const file=e.target.files[0];if(!file)return;const raw=JSON.parse(await file.text()),config=raw.config||raw,validation=await api("/api/validate",config);if(!validation.valid)throw new Error((validation.errors||["Invalid configuration file."]).join("\n"));S.config=normalizeConfig(config);S.presetId="";$("#preset-select").value="";markDirty();renderEditor();toast("Configuration imported. Apply & restart when ready.");e.target.value="";}));
  $("#export-btn").addEventListener("click",()=>{if(S.config)download(S.config,fileName("configuration"));});
  $("#studies-btn").addEventListener("click",()=>navigateEditor("studies"));
  $("#cluster-btn").addEventListener("click",()=>{$("#editor-nav [data-tab='cluster']").click();$("#mobile-workspace-nav [data-mobile-tab='editor']").click();});
  $("#export-run-btn").addEventListener("click",safe(async()=>{await settleRun();download(await api("/api/export"),fileName("run"));toast("Run data exported.");}));
  $("#home-camera").addEventListener("click",()=>{S.camera={yaw:-.42,pitch:.22,zoom:1};queueRender();});
  $("#zoom-in").addEventListener("click",()=>{S.camera.zoom=Math.min(4,S.camera.zoom*1.18);queueRender();});
  $("#zoom-out").addEventListener("click",()=>{S.camera.zoom=Math.max(.35,S.camera.zoom/1.18);queueRender();});
  $("#full-scene").addEventListener("click",safe(async()=>{const el=$("#scene-container");if(document.fullscreenElement)await document.exitFullscreen();else if(el.requestFullscreen)await el.requestFullscreen();else el.classList.toggle("expanded");queueRender();}));
  [["#show-bonds","bonds"],["#show-proteins","proteins"],["#show-labels","labels"],["#show-rna","rna"],["#show-occupancy","occupancy"],["#show-nucleus","nucleus"],["#fit-nucleus","fitNucleus"]].forEach(([id,key])=>$(id).addEventListener("change",e=>{S.view[key]=e.target.checked;queueRender();}));
  let pointer=null;
  $("#scene").addEventListener("pointerdown",e=>{pointer={x:e.clientX,y:e.clientY,total:0};e.currentTarget.setPointerCapture(e.pointerId);e.currentTarget.classList.add("dragging");});
  $("#scene").addEventListener("pointermove",e=>{if(!pointer)return;const dx=e.clientX-pointer.x,dy=e.clientY-pointer.y;pointer.total+=Math.abs(dx)+Math.abs(dy);pointer.x=e.clientX;pointer.y=e.clientY;S.camera.yaw+=dx*.006;S.camera.pitch=Math.max(-1.45,Math.min(1.45,S.camera.pitch+dy*.006));queueRender();});
  $("#scene").addEventListener("pointerup",e=>{if(pointer&&pointer.total<5){const rect=e.currentTarget.getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top;const hits=S.hitTargets.filter(t=>Math.hypot(t.x-x,t.y-y)<t.r).sort((a,b)=>{if(a.type==="body"&&b.type!=="body")return 1;if(b.type==="body"&&a.type!=="body")return-1;if(Math.abs(a.z-b.z)<1e-6){if(a.type==="transcript"&&b.type!=="transcript")return-1;if(b.type==="transcript"&&a.type!=="transcript")return 1;}return b.z-a.z;});selectScene(hits[0]||null,true);}pointer=null;e.currentTarget.classList.remove("dragging");});
  $("#scene").addEventListener("pointercancel",()=>{pointer=null;$("#scene").classList.remove("dragging");});
  $("#scene").addEventListener("wheel",e=>{e.preventDefault();S.camera.zoom=Math.max(.35,Math.min(4,S.camera.zoom*Math.exp(-e.deltaY*.001)));queueRender();},{passive:false});
  $("#scene").addEventListener("keydown",e=>{const keys={ArrowLeft:[-.1,0],ArrowRight:[.1,0],ArrowUp:[0,-.1],ArrowDown:[0,.1]};if(keys[e.key]){e.preventDefault();S.camera.yaw+=keys[e.key][0];S.camera.pitch=Math.max(-1.45,Math.min(1.45,S.camera.pitch+keys[e.key][1]));queueRender();}if(e.key==="Escape")selectScene(null);});
  window.addEventListener("resize",()=>{S.chartsDirty=true;queueRender();});document.addEventListener("fullscreenchange",queueRender);
  new ResizeObserver(()=>{S.chartsDirty=true;queueRender();}).observe($("#scene-container"));
  function canonicalJSON(value){
    if(Array.isArray(value))return `[${value.map(canonicalJSON).join(",")}]`;
    if(value&&typeof value==="object")return `{${Object.keys(value).sort().map(key=>`${JSON.stringify(key)}:${canonicalJSON(value[key])}`).join(",")}}`;
    return JSON.stringify(value);
  }
  async function init(){
    updateControls();
    // Opening or reloading a viewer must never restart the server-side experiment.
    const results=await Promise.allSettled([api("/api/config"),api("/api/presets"),api("/api/capabilities"),api("/api/state")]);
    if(results[0].status!=="fulfilled"||results[3].status!=="fulfilled"){
      const cause=results[0].status!=="fulfilled"?results[0].reason:results[3].reason;
      throw new Error(`Could not restore the current session: ${cause.message||cause}. Check the local server and reload this page. No restart was requested.`);
    }
    let config=results[0].value,snapshot=results[3].value,exported=null;
    // Export is a read-only, coherent view of configuration, state, and recorded history.
    try{exported=await api("/api/export");if(exported.config&&exported.state){config=exported.config;snapshot=exported.state;}}catch(_){/* Current state is usable even when history export is unavailable. */}
    if(!snapshot||!Array.isArray(snapshot.chromatin))throw new Error("The local server has no readable simulation state. Check the server and reload; this viewer has not restarted the experiment.");
    S.config=normalizeConfig(clone(config));S.applied=clone(S.config);
    S.backend=exported?.metadata?.backend||snapshot.engine||"numpy";S.appliedBackend=S.backend;
    S.presetId="";S.appliedPresetId="";S.dirty=false;S.ready=true;
    if(results[1].status==="fulfilled")S.presets=Array.isArray(results[1].value)?results[1].value:results[1].value.presets||[];
    if(results[2].status==="fulfilled")S.capabilities=results[2].value;
    S.validationWarnings=S.capabilities?.limitations||[];
    if(Array.isArray(exported?.history))S.history=exported.history.slice(-1000).map(row=>{
      const rna=(row.genes||[]).reduce((a,g)=>{a.nascent+=Number(g.nascent)||0;a.unprocessed+=Number(g.unprocessed)||0;a.mature+=Number(g.mature)||0;return a;},{nascent:0,unprocessed:0,mature:0});
      return {time:row.time,step:row.step,...(row.metrics||{}),...rna};
    });
    $("#preset-select").innerHTML='<option value="">Custom experiment</option>'+S.presets.map(p=>`<option value="${esc(p.id)}">${esc(p.name)}</option>`).join("");
    if(S.capabilities?.backends?.length)$("#engine-select").innerHTML=S.capabilities.backends.map(b=>`<option value="${esc(b.id)}" ${b.available?"":"disabled"}>${b.id==="numpy"?"Brownian reference":b.id==="polychrom"?"Polychrom / OpenMM":esc(b.id)}${b.available?"":" · unavailable"}</option>`).join("");
    if(!Array.from($("#engine-select").options).some(o=>o.value===S.backend))$("#engine-select").add(new Option(S.backend,S.backend));
    $("#engine-select").value=S.backend;
    $("#experiment-title").textContent=S.applied.name||"Untitled experiment";
    $("#connection-status").innerHTML='<i></i> Local session restored';
    acceptSnapshot(snapshot);renderEditor();renderLegend();updateParticleCount();updateControls();
    // Names are presentation text, not a reliable preset identity. Compare exact model content.
    const appliedSignature=canonicalJSON(S.applied);
    const presetConfigs=await Promise.allSettled(S.presets.map(p=>api(`/api/presets/${encodeURIComponent(p.id)}`)));
    const matchingIndex=presetConfigs.findIndex(result=>result.status==="fulfilled"&&canonicalJSON(normalizeConfig(result.value))===appliedSignature);
    if(!S.dirty&&canonicalJSON(S.applied)===appliedSignature&&matchingIndex>=0){S.presetId=S.presets[matchingIndex].id;S.appliedPresetId=S.presetId;$("#preset-select").value=S.presetId;if(S.tab==="system")renderEditor();}
  }
  safe(init)().catch(()=>{});
})();
