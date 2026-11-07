"use strict";
const test=require("node:test");
const assert=require("node:assert/strict");
const {pairKey,pairCatalog,filteredPairs,contextCount,speciesLayout}=require("../static/interaction-view.js");

test("reversed species names join the configured rule to live molecule counts",()=>{
  const rule={a:"writer",b:"reader",strength:0,on_rate:2};
  const live={a:"reader",b:"writer",molecule_pairs:{total:2},links:{total:5}};
  const pairs=pairCatalog({interactions:[rule]},{pairs:[live]});
  assert.equal(pairs.length,1);
  assert.equal(pairs[0].rule,rule);
  assert.equal(contextCount(pairs[0]),2);
  assert.equal(contextCount(pairs[0],"total","links"),5);
});

test("same-species interactions remain distinct from mixed pairs and identifiers cannot collide",()=>{
  assert.equal(pairKey("A","B"),pairKey("B","A"));
  assert.notEqual(pairKey("A","A"),pairKey("A","B"));
  assert.notEqual(pairKey("A:B","C"),pairKey("A","B:C"));
  const pairs=pairCatalog({interactions:[{a:"A",b:"A"},{a:"A",b:"B"}]},{pairs:[]});
  assert.equal(pairs.length,2);
});

test("missing live data stays unavailable while an observed zero remains zero",()=>{
  const pairs=pairCatalog({interactions:[{a:"A",b:"B"}]},null);
  assert.equal(contextCount(pairs[0],"bound_free"),null);
  const [observed]=pairCatalog({interactions:[]},{pairs:[{a:"A",b:"B",molecule_pairs:{bound_free:0}}]});
  assert.equal(contextCount(observed,"bound_free"),0);
  assert.equal(observed.rule,null);
});

test("species filter retains self pairs and either order without inventing rules",()=>{
  const pairs=pairCatalog({interactions:[{a:"A",b:"A"},{a:"B",b:"A"},{a:"B",b:"C"}]},null);
  assert.deepEqual(filteredPairs(pairs,"A").map(p=>[p.a,p.b]),[["A","A"],["B","A"]]);
  assert.equal(filteredPairs(pairs,"").length,3);
  assert.equal(filteredPairs(pairs,"absent").length,0);
});

test("32 species use compact numbered nodes with non-overlapping hit circles",()=>{
  const nodes=speciesLayout(32);
  assert.equal(nodes.length,32);
  for(let i=0;i<nodes.length;i++){
    const a=nodes[i];assert.equal(a.compact,true);
    assert.ok(a.x-a.radius>=0&&a.x+a.radius<=660&&a.y-a.radius>=0&&a.y+a.radius<=390);
    for(const b of nodes.slice(i+1))assert.ok(Math.hypot(a.x-b.x,a.y-b.y)>a.radius+b.radius+2);
  }
  assert.equal(speciesLayout(6).every(node=>!node.compact&&node.radius===17),true);
});
