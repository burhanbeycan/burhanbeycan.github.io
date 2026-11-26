// Check the literature row matching and the same-concentration crossover.
'use strict';
const assert=require('node:assert/strict');
const path=require('node:path'), fs=require('node:fs');
const root=path.dirname(__dirname);
const web=fs.existsSync(path.join(root,'dist'))?path.join(root,'dist'):root;
global.window=globalThis;
require(path.join(web,'chain-data.js'));require(path.join(web,'chain-pair-data.js'));
const {overlap}=require(path.join(web,'chain-engine.js'));
const paired=new Map(PS_BENZENE_25C.map(row=>[row.id,row]));
assert.equal(paired.size,6);
for(const row of CHAIN_DATA.rows){
  const b=paired.get(row.id);assert.ok(b && b.Rg_nm>row.Rg_nm);
  const ch=overlap(row.Mw_g_mol,row.Rg_nm,1),bz=overlap(b.Mw_g_mol,b.Rg_nm,1);
  const expected=(row.Mw_g_mol/b.Mw_g_mol)*(b.Rg_nm/row.Rg_nm)**3;
  assert.ok(Math.abs(ch.cstar_g_L/bz.cstar_g_L-expected)<1e-12);
  if(row.id==='IK1500-1'){
    assert.ok(Math.abs(ch.cstar_g_L-1.862993)<1e-6 && Math.abs(bz.cstar_g_L-.173803)<1e-6);
    assert.ok(ch.ratio<1 && bz.ratio>1);
  }
}
console.log('PASS paired PS fractions, radius and mass conventions, same-concentration reference example.');
