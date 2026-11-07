/* Applied interaction rules and live molecule-level connections. */
(function(global){
  "use strict";
  const pairKey=(a,b)=>JSON.stringify([String(a),String(b)].sort());
  function pairCatalog(config,state){
    const byKey=new Map();
    for(const rule of config?.interactions||[])byKey.set(pairKey(rule.a,rule.b),{a:rule.a,b:rule.b,rule,live:null});
    for(const pair of state?.pairs||[]){const key=pairKey(pair.a,pair.b),entry=byKey.get(key)||{a:pair.a,b:pair.b,rule:null};entry.live=pair;byKey.set(key,entry);}
    return Array.from(byKey.entries()).map(([key,pair])=>({key,...pair}));
  }
  function filteredPairs(pairs,species){return species?pairs.filter(p=>p.a===species||p.b===species):pairs;}
  function contextCount(pair,context="total",unit="molecule_pairs"){const value=pair?.live?.[unit]?.[context];return Number.isFinite(value)?value:null;}
  function speciesLayout(count){const compact=count>12;return Array.from({length:count},(_,i)=>{const angle=2*Math.PI*i/Math.max(count,1)-Math.PI/2;return {x:count===1?330:330+238*Math.cos(angle),y:count===1?185:195+120*Math.sin(angle),radius:compact?10:17,compact};});}
  const helpers={pairKey,pairCatalog,filteredPairs,contextCount,speciesLayout};
  if(typeof module!=="undefined"&&module.exports)module.exports=helpers;
  if(typeof document==="undefined")return;
  const esc=s=>String(s??"").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const num=(n,d=2)=>n!==null&&n!==undefined&&Number.isFinite(Number(n))?Number(n).toLocaleString(undefined,{maximumFractionDigits:d}):"—";
  const scientific=n=>n!==null&&n!==undefined&&Number.isFinite(Number(n))?Number(n).toExponential(3):"—";
  const measure=n=>n!==null&&n!==undefined&&Number.isFinite(Number(n))?(Number(n)!==0&&Math.abs(Number(n))<0.001?scientific(n):num(n,5)):"—";
  const contexts=[["free_free","Neither directly DNA-bound"],["bound_free","One directly DNA-bound"],["bound_bound","Both directly DNA-bound"]];
  class InteractionView{
    constructor(container,{estimate=async()=>{throw new Error("Physical estimator unavailable.");}}={}){
      this.root=container;this.estimate=estimate;this.config=null;this.state=null;this.pairs=[];this.selected="";this.filter="";this.context="total";this.component="";this.signature="";this.dirty=false;this.estimateGeneration=0;
      container.innerHTML=`<div class="interaction-heading"><div><span class="eyebrow">PROTEIN INTERACTIONS</span><h3>Connections in context</h3><p>Applied rules and the molecular links present in this frame</p></div><span data-role="time" class="count-tag">— t</span></div>
        <div data-role="applied-note" class="interaction-applied-note"></div>
        <div data-role="summary" class="interaction-summary"></div>
        <div class="interaction-controls"><label>Focus on species<select data-role="species-filter" aria-label="Filter interactions by protein species"><option value="">All proteins</option></select></label><label>Live connections shown<select data-role="context-filter" aria-label="Filter interactions by DNA binding context"><option value="total">All binding contexts</option>${contexts.map(([key,title])=>`<option value="${key}">${title}</option>`).join("")}</select></label></div>
        <div class="interaction-network-wrap"><svg data-role="network" class="interaction-network" viewBox="0 0 660 390" role="img" aria-label="Protein species interaction network"></svg></div>
        <section data-role="species-key" class="interaction-species-key" aria-label="Protein species key" hidden><label>Find a species in the key<input data-role="species-search" type="search" placeholder="Search protein name or identifier" aria-label="Search protein species key"></label><p>Numbered nodes keep large networks readable. Select a species below to focus its interactions.</p><div data-role="species-key-list" class="interaction-species-key-list"></div><span data-role="species-key-count" class="small-muted" aria-live="polite"></span></section>
        <div class="interaction-legend"><span><i class="rule-line"></i>Dashed: configured link rule</span><span><i class="live-line"></i>Solid: live molecule pairs · count on line</span><span>Node size is fixed; layout is schematic</span></div>
        <div class="interaction-pair-picker"><label>Inspect protein pair<select data-role="pair" aria-label="Inspect protein interaction pair"></select></label><span data-role="rule-status" class="count-tag"></span></div>
        <div data-role="pair-details"></div>
        <section class="interaction-assemblies"><div class="interaction-section-title"><h4>Live molecular assembly</h4><span>Actual molecules, including higher-order assemblies</span></div><label class="interaction-assembly-picker">Assembly<select data-role="component" aria-label="Inspect live molecular assembly"></select></label><div data-role="assembly"></div></section>
        <section class="interaction-physical"><div class="interaction-section-title"><h4>Physical binding estimate</h4><span>A proposed target, with explicit unit calibration</span></div><p>Assume diffusion-limited association for the selected pair in dilute solution. Supply the physical length and time represented by one simulation unit.</p>
          <button class="text-button" type="button" data-role="use-mapping" hidden>Use the applied model’s fitted length and time ↗</button><form data-role="kinetics-form" class="interaction-kinetics-form"><label class="interaction-kinetics-target">Binding target<select data-role="kinetics-target" name="binding_target" aria-label="Physical binding estimate target"><option value="pair">Selected protein pair</option></select></label><label>Length per reduced unit · nm<input name="length_nm" type="number" min="0.000000001" step="any" required placeholder="Your length calibration" aria-label="Physical length per reduced unit in nanometers"></label><label>Time per reduced unit · s<input name="time_s" type="number" min="0.000000001" step="any" required placeholder="Your time calibration" aria-label="Physical time per reduced unit in seconds"></label><label>Target Kd · μM<input name="target_kd_um" type="number" min="0.000000001" step="any" required value="1" aria-label="Target dissociation constant in micromolar"></label><label>Capture efficiency · 0 to 1<input name="efficiency" type="number" min="0.000000001" max="1" step="any" required value="1" aria-label="Diffusion capture efficiency"></label><button class="secondary-button" type="submit">Estimate physical rates</button></form>
          <p class="interaction-physical-note">Efficiency 1 is the ideal absorbing-sphere limit. This estimate does not change the simulation. The current encounter hazard and spring strength are not measured kon or Kd.</p><div data-role="kinetics-result" aria-live="polite"></div>
        </section>`;
      this.find=role=>this.root.querySelector(`[data-role="${role}"]`);
      this.find("species-filter").addEventListener("change",e=>{this.filter=e.target.value;this.renderPairPicker();this.renderNetwork();this.renderPair();});
      this.find("context-filter").addEventListener("change",e=>{this.context=e.target.value;this.renderNetwork();});
      this.find("species-search").addEventListener("input",()=>this.renderSpeciesKey());
      this.find("species-key-list").addEventListener("click",e=>{const button=e.target.closest("[data-key-species]");if(button)this.focusSpecies(button.dataset.keySpecies);});
      this.find("pair").addEventListener("change",e=>{this.selected=e.target.value;this.clearEstimate();this.renderNetwork();this.renderPair();});
      this.find("component").addEventListener("change",e=>{this.component=e.target.value;this.renderAssembly();});
      this.find("network").addEventListener("click",e=>{const pair=e.target.closest("[data-pair-key]"),species=e.target.closest("[data-species-id]");if(pair){this.selected=pair.dataset.pairKey;this.find("pair").value=this.selected;this.clearEstimate();this.renderNetwork();this.renderPair();}else if(species){this.filter=this.filter===species.dataset.speciesId?"":species.dataset.speciesId;this.find("species-filter").value=this.filter;this.renderPairPicker();this.renderNetwork();this.renderPair();}});
      this.find("network").addEventListener("keydown",e=>{if(e.key==="Enter"||e.key===" "){e.preventDefault();e.target.dispatchEvent(new MouseEvent("click",{bubbles:true}));}});
      this.find("use-mapping").addEventListener("click",()=>{const p=this.config?.calibration;if(p?.status!=="fitted")return;const form=this.find("kinetics-form");form.elements.length_nm.value=p.length_nm;form.elements.time_s.value=p.time_s;this.clearEstimate();});
      this.find("kinetics-form").addEventListener("input",()=>this.clearEstimate());
      this.find("kinetics-target").addEventListener("change",()=>{this.clearEstimate();this.updateEstimateAvailability();});
      this.find("kinetics-form").addEventListener("submit",e=>{e.preventDefault();this.calculate();});
    }
    update(snapshot,config,{dirty=false}={}){
      this.find("use-mapping").hidden=config?.calibration?.status!=="fitted";
      this.snapshot=snapshot;this.config=config;this.state=snapshot?.interaction_state||null;this.dirty=dirty;
      this.pairs=pairCatalog(config,this.state);
      const signature=JSON.stringify([config?.proteins,config?.interactions,config?.chromatin]);
      if(signature!==this.signature){this.signature=signature;this.clearEstimate();const picker=this.find("species-filter");picker.innerHTML='<option value="">All proteins</option>'+(config?.proteins||[]).map(p=>`<option value="${esc(p.id)}">${esc(p.name)}</option>`).join("");if(!(config?.proteins||[]).some(p=>p.id===this.filter))this.filter="";picker.value=this.filter;const target=this.find("kinetics-target"),value=target.value;target.innerHTML='<option value="pair">Selected protein pair · dilute solution</option>'+(config?.proteins||[]).map(p=>`<option value="dna:${esc(p.id)}">${esc(p.name)} → chromatin · isolated bead proxy</option>`).join("");target.value=Array.from(target.options).some(o=>o.value===value)?value:"pair";}
      this.setDirty(dirty);this.find("time").textContent=`${num(snapshot?.time,3)} t · step ${num(snapshot?.step,0)}`;
      const totals=this.state?.totals;
      this.find("summary").innerHTML=totals?[[totals.direct_bound,"directly DNA-bound"],[totals.indirect_bound,"indirectly tethered"],[totals.free,"not DNA-tethered"],[totals.assemblies,"assemblies ≥2 molecules"]].map(([value,label])=>`<div><strong>${num(value,0)}</strong><span>${label}</span></div>`).join(""):'<p class="interaction-empty">Live connection data is unavailable from this server. Applied rules remain visible.</p>';
      this.renderPairPicker();this.renderNetwork();this.renderPair();this.renderAssemblies();
    }
    setDirty(dirty){this.dirty=dirty;this.find("applied-note").textContent=dirty?"Showing the applied model. Draft edits appear after Apply & restart.":"Showing the applied model. Counts update with each simulation frame.";}
    focusSpecies(id){if(!(this.config?.proteins||[]).some(p=>p.id===id))return;this.filter=id;this.find("species-filter").value=id;this.renderPairPicker();this.renderNetwork();this.renderPair();}
    species(id){return (this.config?.proteins||[]).find(p=>p.id===id)||{id,name:id,color:"#819096"};}
    selectedPair(){return this.pairs.find(p=>p.key===this.selected);}
    pairName(pair){return `${this.species(pair.a).name} ↔ ${this.species(pair.b).name}`;}
    renderPairPicker(){
      const visible=filteredPairs(this.pairs,this.filter),old=this.selected;
      if(!visible.some(p=>p.key===this.selected))this.selected=visible[0]?.key||"";
      const signature=JSON.stringify(visible.map(p=>[p.key,this.pairName(p)]));
      if(signature!==this.pairOptionsSignature){this.pairOptionsSignature=signature;this.find("pair").innerHTML=visible.length?visible.map(p=>`<option value="${esc(p.key)}">${esc(this.pairName(p))}</option>`).join(""):'<option value="">No interaction rules for this selection</option>';}
      this.find("pair").value=this.selected;this.find("pair").disabled=!visible.length;
      if(old!==this.selected)this.clearEstimate();
    }
    renderNetwork(){
      const all=this.config?.proteins||[],pairs=filteredPairs(this.pairs,this.filter),ids=new Set(pairs.flatMap(p=>[p.a,p.b]));if(this.filter)ids.add(this.filter);
      const species=this.filter?all.filter(p=>ids.has(p.id)):all,positions=new Map(),layout=speciesLayout(species.length),numbers=new Map(all.map((p,i)=>[p.id,i+1]));
      species.forEach((p,i)=>positions.set(p.id,layout[i]));
      let svg='<title>Applied protein interaction rules and current molecular connections</title><desc>Dashed lines indicate configured possible interactions. Solid lines and numbers indicate distinct molecule pairs linked in the current frame. Select a line to inspect its pair.</desc>';
      for(const pair of pairs){const a=positions.get(pair.a),b=positions.get(pair.b);if(!a||!b)continue;const count=contextCount(pair,this.context),active=count!==null&&count>0,configured=!!pair.rule,selected=pair.key===this.selected;if(!configured&&!active)continue;const self=pair.a===pair.b,path=self?`M ${a.x-9} ${a.y-15} C ${a.x-65} ${a.y-75}, ${a.x+65} ${a.y-75}, ${a.x+9} ${a.y-15}`:`M ${a.x} ${a.y} L ${b.x} ${b.y}`,tx=self?a.x:(a.x+b.x)/2,ty=self?a.y-56:(a.y+b.y)/2;
        svg+=`<g data-pair-key="${esc(pair.key)}" tabindex="0" role="button" aria-label="Inspect ${esc(this.pairName(pair))}: ${num(count,0)} live molecule pairs" class="interaction-network-edge"><title>${esc(this.pairName(pair))} · ${num(count,0)} live pairs${configured?" · configured":""}</title><path d="${path}" fill="none" stroke="transparent" stroke-width="17"/><path d="${path}" fill="none" stroke="${selected?"#b98549":active?"#197f79":"#b5c3c0"}" stroke-width="${selected?3:active?2.5:1.5}" ${active?'':'stroke-dasharray="5 5"'}/>${active?`<rect x="${tx-12}" y="${ty-9}" width="24" height="18" rx="7" fill="white"/><text x="${tx}" y="${ty+4}" text-anchor="middle" fill="#197f79" font-size="11" font-weight="600">${num(count,0)}</text>`:""}</g>`;
      }
      for(const p of species){const {x,y,radius,compact}=positions.get(p.id),s=this.state?.species?.find(s=>s.id===p.id),number=numbers.get(p.id);svg+=`<g data-species-id="${esc(p.id)}" tabindex="0" role="button" aria-label="Focus interactions on ${esc(p.name)}" class="interaction-network-node"><title>${compact?`${number} · `:""}${esc(p.name)} · ${num(s?.total??p.count,0)} molecules · ${num(s?.direct_bound,0)} directly bound · ${num(s?.indirect_bound,0)} indirectly tethered</title><circle cx="${x}" cy="${y}" r="${radius}" fill="${esc(p.color)}" stroke="${this.filter===p.id?"#283a46":"white"}" stroke-width="${compact?2:3}"/><text x="${x}" y="${compact?y+3.5:y+33}" text-anchor="middle" fill="${compact?"#fff":"#41594e"}" font-size="${compact?10:11}" ${compact?'font-weight="650"':''}>${compact?number:esc(p.name.length>21?p.name.slice(0,19)+"…":p.name)}</text></g>`;}
      if(!species.length)svg+='<text x="330" y="180" text-anchor="middle" fill="#819096">No protein species in the applied model</text>';
      this.find("network").innerHTML=svg;
      this.renderSpeciesKey();
    }
    renderSpeciesKey(){
      const all=this.config?.proteins||[],container=this.find("species-key");container.hidden=all.length<=12;if(container.hidden)return;
      const query=this.find("species-search").value.trim().toLocaleLowerCase(),visible=all.map((p,i)=>({protein:p,number:i+1})).filter(({protein:p,number})=>!query||`${p.name} ${p.id} ${number}`.toLocaleLowerCase().includes(query));
      const signature=JSON.stringify([visible.map(({protein:p,number})=>[p.id,p.name,p.color,number]),this.filter]);
      if(signature!==this.speciesKeySignature){this.speciesKeySignature=signature;const focused=document.activeElement?.dataset?.keySpecies,list=this.find("species-key-list");list.innerHTML=visible.length?visible.map(({protein:p,number})=>`<button type="button" data-key-species="${esc(p.id)}" aria-pressed="${this.filter===p.id}" aria-label="Focus interactions on ${esc(p.name)}"><span style="background:${esc(p.color)}">${number}</span><strong>${esc(p.name)}</strong></button>`).join(""):'<p class="interaction-key-empty">No matching species. Try another name or identifier.</p>';if(focused)Array.from(list.querySelectorAll("button")).find(button=>button.dataset.keySpecies===focused)?.focus();}
      this.find("species-key-count").textContent=`${visible.length} of ${all.length} species`;
    }
    contextDiagram(key,a,b){
      const boundA=key!=="free_free",boundB=key==="bound_bound";
      return `<svg viewBox="0 0 220 120" role="img" aria-label="Schematic: ${esc(contexts.find(c=>c[0]===key)[1])}"><path d="M20 98 Q50 85 80 98 T140 98 T200 98" fill="none" stroke="#859daf" stroke-width="3"/><text x="110" y="115" text-anchor="middle" fill="#6a8092" font-size="9">chromatin / DNA</text>${boundA?'<path d="M70 52 L70 94" stroke="#718fa2" stroke-width="3"/>':""}${boundB?'<path d="M150 52 L150 94" stroke="#718fa2" stroke-width="3"/>':""}<path d="M83 40 L137 40" stroke="#197f79" stroke-width="3"/><circle cx="70" cy="40" r="17" fill="${esc(a.color)}"/><circle cx="150" cy="40" r="17" fill="${esc(b.color)}"/><text x="70" y="44" text-anchor="middle" font-size="11" fill="white">A</text><text x="150" y="44" text-anchor="middle" font-size="11" fill="white">B</text></svg>`;
    }
    renderPair(){
      const pair=this.selectedPair();this.find("rule-status").textContent=pair?(pair.rule?Number(pair.rule.on_rate)>0?"APPLIED RULE":"ASSOCIATION OFF":"LIVE CONNECTION"):"NO PAIR";
      this.updateEstimateAvailability();
      if(!pair){this.find("pair-details").innerHTML='<p class="interaction-empty">Add interaction rules in the model builder, then apply the model to inspect protein pairs here.</p>';return;}
      const a=this.species(pair.a),b=this.species(pair.b),rule=pair.rule;
      this.find("pair-details").innerHTML=`<div class="interaction-pair-title"><span><i style="background:${esc(a.color)}"></i>A · ${esc(a.name)}</span><span><i style="background:${esc(b.color)}"></i>B · ${esc(b.name)}${pair.a===pair.b?" · distinct molecules of the same species":""}</span></div><div class="interaction-context-grid">${contexts.map(([key,title])=>`<article class="interaction-context-card"><h4>${title}</h4>${this.contextDiagram(key,a,b)}<strong>${num(contextCount(pair,key),0)} <small>molecule pairs</small></strong><span>${num(contextCount(pair,key,"links"),0)} individual links</span></article>`).join("")}</div><p class="interaction-definition">These three contexts count <strong>direct DNA attachment</strong> of the two interacting molecules. In the middle diagram, either partner may be attached; one example is shown. A molecule without a direct attachment can still be tethered through other proteins. Multiple IDR domain links between the same two molecules count as one molecule pair.</p>${rule?`<div class="interaction-rule-values">${[["Link spring strength",num(rule.strength,3)],["Encounter radius",`${num(rule.radius,3)} length`],["Encounter hazard",`${measure(rule.on_rate)} / time`],["Link off-rate",`${measure(rule.off_rate)} / time`]].map(([label,value])=>`<div><span>${label}</span><strong>${value}</strong></div>`).join("")}</div><p class="interaction-definition">The same applied pair rule is used in all three contexts. Binding changes motion through mechanical links; this view does not assign new context-specific affinities. Units above are reduced units.</p>`:""}`;
    }
    renderAssemblies(){
      const components=(this.state?.components||[]).filter(c=>c.size>1).sort((a,b)=>b.size-a.size||a.id-b.id),picker=this.find("component");
      if(!components.some(c=>String(c.id)===this.component))this.component=components.length?String(components[0].id):"";
      picker.innerHTML=components.length?components.map(c=>`<option value="${c.id}">Assembly ${c.id} · ${c.size} molecules · ${c.chromatin_bound?"DNA-tethered":"not DNA-tethered"}</option>`).join(""):'<option value="">No protein assembly in this frame</option>';picker.value=this.component;picker.disabled=!components.length;this.renderAssembly();
    }
    renderAssembly(){
      const c=this.state?.components?.find(c=>String(c.id)===this.component),container=this.find("assembly");
      if(!c||c.size<2){container.innerHTML=`<p class="interaction-empty">${this.state?"No connected multi-protein assembly is present in this frame. Pair diagrams above illustrate the binding contexts even when their current counts are zero.":"Live assembly data is unavailable."}</p>`;return;}
      const ids=c.molecule_ids.slice(0,48),idSet=new Set(ids),molecules=new Map((this.state.molecules||[]).map(m=>[m.id,m])),positions=new Map();ids.forEach((id,i)=>{const angle=i/ids.length*2*Math.PI-Math.PI/2;positions.set(id,{x:300+205*Math.cos(angle),y:157+108*Math.sin(angle)});});
      let svg='<title>Live protein assembly. Each node is one molecule.</title><path d="M30 308 Q90 298 150 308 T270 308 T390 308 T570 308" stroke="#859daf" stroke-width="3" fill="none"/><text x="300" y="330" text-anchor="middle" font-size="10" fill="#718fa2">chromatin · direct attachments shown in blue</text>';
      for(const id of ids){const m=molecules.get(id),p=positions.get(id);if(m?.binding_status==="direct")svg+=`<path d="M${p.x} ${p.y} L${p.x} 308" stroke="#a4b9c8" stroke-width="1.5"/>`;}
      for(const link of this.state.molecule_links||[]){if(!idSet.has(link.a)||!idSet.has(link.b))continue;const a=positions.get(link.a),b=positions.get(link.b);svg+=`<path d="M${a.x} ${a.y} L${b.x} ${b.y}" stroke="#197f79" stroke-width="${Math.min(5,1+link.links)}"><title>${link.a} ↔ ${link.b}: ${link.links} individual links</title></path>`;}
      for(const id of ids){const m=molecules.get(id),p=positions.get(id),species=this.species(m?.species);svg+=`<g><title>Molecule ${id} · ${esc(species.name)} · ${m?.binding_status==="direct"?"directly DNA-bound":m?.binding_status==="indirect"?"indirectly DNA-tethered":"not DNA-tethered"}</title>${m?.binding_status==="indirect"?`<circle cx="${p.x}" cy="${p.y}" r="18" fill="none" stroke="#b98549" stroke-width="2" stroke-dasharray="3 3"/>`:""}<circle cx="${p.x}" cy="${p.y}" r="13" fill="${esc(species.color)}" stroke="white" stroke-width="2"/><text x="${p.x}" y="${p.y+4}" text-anchor="middle" font-size="9" fill="white">${id}</text></g>`;}
      container.innerHTML=`<svg class="interaction-assembly-network" viewBox="0 0 600 340" role="img" aria-label="Current assembly ${c.id} molecular links">${svg}</svg><div class="interaction-assembly-stats"><span><strong>${num(c.size,0)}</strong> molecules</span><span><strong>${num(c.direct_bound,0)}</strong> direct</span><span><strong>${num(c.indirect_bound,0)}</strong> indirect</span><span><strong>${num(c.protein_links,0)}</strong> links</span></div><p class="interaction-definition">Node labels are molecule IDs. Gold dashed rings mark indirect tethering. Layout is schematic; protein connections are current.${c.size>48?` Displaying the first 48 of ${c.size} molecules; totals include the complete assembly.`:""} Assembly IDs describe this frame, not persistent lineages.</p>`;
    }
    clearEstimate(){this.estimateGeneration++;if(this.find)this.find("kinetics-result").innerHTML="";}
    estimateTarget(){const value=this.find("kinetics-target").value;if(value.startsWith("dna:"))return {a:value.slice(4),b:"chromatin"};const pair=this.selectedPair();return pair?.rule?{a:pair.a,b:pair.b}:null;}
    updateEstimateAvailability(){this.find("kinetics-form").querySelector('button').disabled=!this.estimateTarget();}
    async calculate(){
      const target=this.estimateTarget();if(!target)return;
      const form=this.find("kinetics-form"),values=Object.fromEntries(new FormData(form)),request={...target,...Object.fromEntries(Object.entries(values).filter(([key])=>key!=="binding_target").map(([key,value])=>[key,Number(value)]))};
      const generation=++this.estimateGeneration,output=this.find("kinetics-result");output.innerHTML='<p class="interaction-definition">Calculating physical rates…</p>';
      try{const result=await this.estimate(request);if(generation!==this.estimateGeneration)return;
        output.innerHTML=`<div class="interaction-physical-results"><div><span>Pair encounter limit · diffusion ceiling</span><strong>${scientific(result.kon_limit_M_inv_s)} <small>M⁻¹ s⁻¹</small></strong></div><div><span>Pair kernel · chosen efficiency</span><strong>${scientific(result.kon_M_inv_s)} <small>M⁻¹ s⁻¹</small></strong></div><div><span>Complex formation coefficient</span><strong>${scientific(result.macroscopic_kon_M_inv_s)} <small>M⁻¹ s⁻¹</small></strong></div><div><span>Nominal koff = Kd × formation coefficient</span><strong>${scientific(result.nominal_koff_s_inv)} <small>s⁻¹</small></strong></div><div><span>Nominal koff · reduced units</span><strong>${scientific(result.nominal_koff_reduced)} <small>/ time</small></strong></div></div><p class="interaction-definition">Free molecule D: A ${measure(result.physical_D_a_um2_s)}, B ${measure(result.physical_D_b_um2_s)} μm²/s. Capture radius: ${measure(result.physical_radius_um*1000)} nm. ${result.pair?.kind==="protein_chromatin_proxy"?"B is an isolated chromatin bead proxy using its bare D and the model capture radius. A tethered polymer is not an independent sphere.":"This uses the applied pair radius and free-molecule diffusion values."}</p><p class="interaction-physical-note">${result.pair?.homodimer?"For A + A ⇌ A₂, Kd = [A]²/[A₂]. The complex formation coefficient is half the pair encounter kernel because identical pairs are unordered.":"For A + B ⇌ AB, Kd = [A][B]/[AB] = koff / complex formation coefficient."} Bound, multivalent, confined, and rebinding systems require calibration before these numbers describe their observed affinity.</p>${Array.isArray(result.limitations)&&result.limitations.length?`<details class="interaction-limitations"><summary>Assumptions and limits of this estimate</summary><ul>${result.limitations.map(item=>`<li>${esc(item)}</li>`).join("")}</ul></details>`:""}`;
      }catch(error){if(generation===this.estimateGeneration)output.innerHTML=`<p class="interaction-calculation-error" role="alert">${esc(error.message||error)}</p>`;}
    }
  }
  global.ChromatinInteractionView=InteractionView;
})(typeof window!=="undefined"?window:globalThis);
