'use strict';
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const base = path.join(__dirname,'..');
const web = fs.existsSync(path.join(base,'dist','index.html')) ? path.join(base,'dist') : base;
const fh = require(path.join(web,'phase-engine.js'));
const osm = require(path.join(web,'osmotic-engine.js'));
const near = (a,b,rel=1e-8,abs=0)=>assert.ok(Math.abs(a-b)<=abs+rel*Math.max(Math.abs(a),Math.abs(b)),a+' ≠ '+b);
const Vs=100, v=.9, Mn=100*Vs/v,T=300;
const e=osm.state(T,Mn,v,Vs,{A:.6,B:0,C:.12});
const zero=osm.evaluate(0,e);
near(zero.pressure,0); near(zero.reduced,osm.R*T/Mn); near(e.N,100);
near(e.A2,.02*v*v/Vs); // Composition slope has positive sign in the mass-basis A2.
const dilute=osm.evaluate(1e-4,e), cml=dilute.c/1000;
const virial=(dilute.pressure/(1000*osm.R*T)-cml/Mn)/cml**2;
near(virial,e.A2,1e-5); near(dilute.phi,9e-8); near(dilute.ideal,osm.R*T*1e-4/Mn);
const at=osm.evaluate(22,e);
near(at.pressure,-1000*osm.R*T/Vs*(fh.f(at.phi,100,e.chi0,e.C)-at.phi*fh.df(at.phi,100,e.chi0,e.C)));
const demix=osm.state(300,Mn,v,Vs,{A:.65,B:0,C:0});
const bins=demix.eq.bin; assert.ok(bins && bins.length===2);
const p1=osm.evaluate(.02*1000/v,demix),p2=osm.evaluate(.16*1000/v,demix);
near(p1.pressure,p2.pressure); assert.equal(p1.phase.count,2);
for(const phi of bins){const endpoint=osm.evaluate(phi*1000/v,demix);near(endpoint.pressure,p1.pressure,1e-7);}
assert.ok(osm.evaluate(.5*1000/v,demix).pressure>p2.pressure);
assert.throws(()=>osm.evaluate(-1,e));assert.throws(()=>osm.evaluate(11000,e));
console.log('PASS osmotic pressure and units, zero/dilute limit, A2 composition derivative, agreement with free energy and pressure plateau.');
// External membrane series: preserve the observed concentration basis and the measured Mn.
global.window={};
require(path.join(web,'osmotic-pressure-data.js'));
require(path.join(web,'transfer-results.js'));
const source=window.OSMOTIC_PRESSURE_DATA;
const group=window.TRANSFER_RESULTS;
assert.equal(source.Mn_g_mol,72000);
assert.equal(source.observations.length,8);
const measurements=source.observations.map(row=>{
  const expected=1e-6*osm.R*row.T_K*row.c_g_L*row.Q1e6;
  near(row.Pi_kPa,expected,0,1e-8);
  return row;
});
const errors={};
for(const model of group.models){
  const prediction=measurements.map(row=>{
    const state=osm.state(row.T_K,source.Mn_g_mol,group.constants.vbar_ml_g,group.constants.V_s_ml_mol,model);
    const res=osm.evaluate(row.c_g_L,state);
    assert.equal(res.phase.count,1);
    return Math.abs(res.pressure-row.Pi_kPa);
  });
  errors[model.id]=prediction.reduce((a,b)=>a+b,0)/prediction.length;
}
near(errors.cloud,7.63631,3e-5);near(errors.combined,10.06294,3e-5);
console.log('PASS eight experimental pressure conversions and full-pressure external errors (kPa)',errors);
