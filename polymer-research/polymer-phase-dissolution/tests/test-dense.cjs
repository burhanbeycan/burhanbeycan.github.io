// node tests/test-dense.cjs [web-directory]
'use strict';
const {strict:assert} = require('node:assert');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const site = path.resolve(__dirname,'..');
const web = path.resolve(process.argv[2] || (fs.existsSync(path.join(site,'dist'))?path.join(site,'dist'):site));
const fh = require(path.join(web,'dense-engine.js'));
const osm = require(path.join(web,'osmotic-engine.js'));
const ctx = {window:{}};
vm.runInNewContext(fs.readFileSync(path.join(web,'dense-results.js'),'utf8'),ctx);
const r = ctx.window.DENSE_RESULTS;
const {Mn_g_mol:Mn,vbar_ml_g:v,Vs_ml_mol:Vs}=r.sample;
const close=(a,b,tol)=>assert.ok(Math.abs(a-b)<=tol,`${a} vs ${b} (tol ${tol})`);
const state=(T,D)=>fh.state(T,Mn,v,Vs,r.base,D);
for(const T of [297.15,307.15,317.15]){
  const orig=osm.state(T,Mn,v,Vs,r.base), s=state(T,0);
  assert.deepEqual(s.eq.bin,orig.eq.bin);
  for(const c of [0,1e-8,1,100,259,600]){
    close(fh.evaluate(c,s).pressure,osm.evaluate(c,orig).pressure,1e-11);
  }
  for(const d of [r.bounds[0],-0.21,-0.08,r.D_guard_min,r.bounds[1]]){
    const a=state(T,d),zero=state(T,0);
    assert.equal(a.A2,zero.A2);
    // In the homogeneous virial expansion with c in g/mL, Π/(1000 RT) =
    // c/Mn + A2*c² + A3*c³ + ...; use the stable solvent-potential expansion.
    const pSmall=5e-5,cMl=pSmall/v;
    const apparentA3=(-fh.solventPotential(pSmall,a)/Vs-cMl/Mn-a.A2*cMl**2)/cMl**3;
    const analyticA3=v**3/Vs*(1/3-2*a.C+2*a.D);
    close(apparentA3,analyticA3,2e-6);
    const tinyC=Math.min(1e-7,(a.eq.bin?.[0] ?? 1)*100/v);
    close(fh.evaluate(tinyC,a).pressure/tinyC,osm.R*T/Mn,1e-8);
    for(const p of [.001,.02,.15,.55]){
      close((fh.energy(p+1e-6,a)-fh.energy(p-1e-6,a))/2e-6,fh.derivative(p,a),2e-9);
      close((fh.derivative(p+1e-6,a)-fh.derivative(p-1e-6,a))/2e-6,fh.curvature(p,a),2e-6);
      close(fh.solventPotential(p,a),fh.energy(p,a)-p*fh.derivative(p,a),1e-13);
    }
    // Independent 2001-point support-line audit across the complete composition range.
    if(a.eq.bin){
      const [l,b]=a.eq.bin;
      close(fh.derivative(l,a),fh.derivative(b,a),1e-10);
      close(fh.solventPotential(l,a),fh.solventPotential(b,a),1e-12);
      const m=fh.derivative(b,a), intercept=fh.solventPotential(b,a);
      for(let i=0;i<=2000;i++){
        const p=i/2000; assert.ok(fh.energy(p,a)-m*p-intercept>=-1e-12);
      }
      const inside=[.1,.5,.9].map(t=>fh.evaluate((l+(b-l)*t)*1000/v,a));
      assert.ok(inside.every(x=>x.phase.count===2 && x.pressure===inside[0].pressure));
    }
  }
}
assert.throws(()=>state(298,fh.maxD+.001),RangeError);
const idsWeExpect=['hold24','hold34','high','hold44'];
assert.deepEqual([...r.folds.map(f=>f.id)],idsWeExpect);
const high=r.folds.find(f=>f.id==='high');
assert.ok(high.rows.filter(row=>row.role==='test').every(row=>!high.rows.some(x=>x.T_K===row.T_K && x.c_g_L>row.c_g_L)));
for(const f of r.folds){
  assert.equal(new Set([...f.train_ids,...f.test_ids]).size,f.train_ids.length+f.test_ids.length);
  assert.ok(f.D>r.bounds[0] && f.D<r.bounds[1]);
  assert.ok(f.A3<0 && f.A3_guard>=0 && r.base.A3>0);
  for (const row of f.rows){
    for(const [key,d] of [['base',0],['extended',f.D],['guard',f.D_guard]]){
      const p=fh.evaluate(row.c_g_L,state(row.T_K,d));
      close(p.pressure,row[key+'_kPa'],1e-10);assert.equal(p.phase.count,row[key+'_phases']);
    }
  }
  const held=f.rows.filter(x=>f.score_ids.includes(x.id));
  for(const key of ['base','extended','guard']){
    const mean=held.reduce((n,row)=>n+Math.abs(row[key+'_kPa']-row.observed_kPa),0)/held.length;
    close(mean,f.metrics[key].mae,1e-10);
  }
}
const ext=r.folds.find(f=>f.id==='hold44');
assert.equal(ext.train_ids.length,8);
assert.equal(ext.test_ids.length,3);
assert.ok(ext.score_ids.length>=1 && ext.score_ids.length<=3);
assert.ok(ext.rows.filter(x=>x.role==='test').every(row=>row.T_K===317.15));
assert.ok(ext.rows.filter(x=>x.role==='train').every(row=>row.T_K<317.15));
const mid=ext.rows.find(x=>x.T_K===317.15 && x.c_g_L===178);
close(mid.observed_kPa,317.15*osm.R*178*28.1/1e6,1e-7);
close(ext.D,-.0743682,2e-5);
console.log('PASS base recovery; derivatives, dilute terms, global tie-line support, plateau, held-out split independence, finite pressure/metrics and 44 °C unit conversions.');
