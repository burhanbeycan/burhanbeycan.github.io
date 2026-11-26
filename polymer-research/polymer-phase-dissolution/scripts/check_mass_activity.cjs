// Score Table III fractions I and III with the frozen pressure fit from fraction II.
// Usage: node scripts/check_mass_activity.cjs [dist-or-flat-web-dir] [input.csv] [output.json]
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');

const baseDir = path.resolve(__dirname, '..');
const web = path.resolve(process.argv[2] || (fs.existsSync(path.join(baseDir, 'dist')) ? path.join(baseDir, 'dist') : baseDir));
const input = process.argv[3] || path.join(web, 'transfer-inputs', 'ps_cyclohexane_mass_activity.csv');
const output = process.argv[4] || path.join(web, 'transfer-inputs', 'mass-activity-results.json');
const table = fs.readFileSync(input, 'utf8').trim().split(/\r?\n/);
assert.equal(table.shift(), 'fraction,Mn_g_mol,T_K,phi,reference_100x_TPM,observed_minus_ln_a');
const rows = table.map(line => {
  const [fraction, ...fields] = line.split(',');
  const [Mn_g_mol,T_K,phi,reference_100x_TPM,observed_minus_ln_a] = fields.map(Number);
  assert.ok(fields.length === 5 && fields.every(x => x.trim() !== '' && Number.isFinite(Number(x))));
  assert.ok((fraction === 'I' && Mn_g_mol === 25100) || (fraction === 'III' && Mn_g_mol === 440000));
  assert.ok([307.15,317.15].includes(T_K) && phi>0 && phi<1 && reference_100x_TPM>0 && observed_minus_ln_a>0);
  // Independent internal transcription check: authors' smooth calibration of the TPM reference.
  // Table values and the printed interpolation both carry rounding.
  const xTPM = reference_100x_TPM/100;
  const [b,h] = T_K<310 ? [3.50,8.8] : [3.18,7.2];
  const minusLogRef = -Math.log((1-xTPM)*(1+b*xTPM*xTPM-h*xTPM**3));
  assert.ok(Math.abs(minusLogRef-observed_minus_ln_a)<0.0006, 'Reference-calibration discrepancy');
  return {fraction,Mn_g_mol,T_K,phi,reference_100x_TPM,observed_minus_ln_a};
});
assert.deepEqual(['I','III'].map(f => rows.filter(r => r.fraction === f).length), [4,9]);
assert.equal(new Set(rows.map(r => [r.fraction,r.T_K,r.phi].join('/'))).size, rows.length);
const ctx = {window:{}};
vm.runInNewContext(fs.readFileSync(path.join(web, 'dense-results.js'), 'utf8'), ctx);
const ref = ctx.window.DENSE_RESULTS;
const fit = ref.folds.find(row => row.id === 'hold44');
assert.equal(fit.train_ids.length,8);
const Dense = require(path.join(web,'dense-engine.js'));
const {vbar_ml_g: v, Vs_ml_mol: Vs} = ref.sample;
const cache = new Map();
const coefficients = {baseline:0,extended:fit.D,guarded:fit.D_guard};
function predict(T, Mn, phi, d) {
  const key = [T,Mn,d].join('/');
  if (!cache.has(key)) cache.set(key, Dense.state(T,Mn,v,Vs,ref.base,d));
  const q = Dense.evaluate(1000*phi/v,cache.get(key));
  const minus_ln_a = q.pressure*Vs/(1000*8.31446261815324*T);
  return {minus_ln_a,a:Math.exp(-minus_ln_a),model_phase_count:q.phase.count};
}
function withPredictions(row, Mn) {
  return {...row,Mn_g_mol:Mn,observed_a:Math.exp(-row.observed_minus_ln_a),
    ...Object.fromEntries(Object.entries(coefficients).map(([id,d])=>[id,predict(row.T_K,Mn,row.phi,d)]))};
}
function score(records) {
  return {n:records.length,...Object.fromEntries(Object.keys(coefficients).map(id=>{
    function mae(field) {return records.reduce((sum,r) => sum+Math.abs(r[id][field]-r['observed_'+field]),0)/records.length;}
    return [id,{mae_minus_ln_a:mae('minus_ln_a'),mae_activity:mae('a')}];
  }))};
}
const scored = rows.map(row => withPredictions(row,row.Mn_g_mol));
const fracI = scored.filter(r => r.fraction === 'I');
const alternativeMass = fracI.map(r => withPredictions(r,25900));
// All recorded compositions lie on homogeneous equilibrium branches for these frozen fits.
for (const row of scored) {
  for (const [id,d] of Object.entries(coefficients)) {
    assert.equal(row[id].model_phase_count, 1);
    const p=row.phi, n=row.Mn_g_mol*v/Vs, chi0=ref.base.A+ref.base.B/row.T_K;
    const homogeneous=-(Math.log1p(-p)+p-p/n+p*p*(chi0-ref.base.C+2*ref.base.C*p)+d*p**3*(3*p-2));
    assert.ok(Math.abs(row[id].minus_ln_a-homogeneous)<1e-12);
  }
}
const result = {
  method:'Same-study transfer from pressure-fit fraction II to additional fractions I and III. No refit. Model phase counts are predictions.',
  source_doi:'10.1021/ja01517a022', table:'III', page:1864, pressure_fit_n:8,
  D:coefficients, scores:{I:score(fracI),III:score(scored.filter(r=>r.fraction==='III')),pooled:score(scored)},
  fraction_I_sensitivity:{alternative_Mn_g_mol:25900,
    max_model_activity_change: Math.max(...alternativeMass.flatMap((r,i)=>Object.keys(coefficients).map(id=>Math.abs(r[id].a-fracI[i][id].a)))),
    scores:score(alternativeMass)},rows:scored
};
fs.writeFileSync(output,JSON.stringify(result,null,2)+'\n');
console.log(JSON.stringify({scores:result.scores,sensitivity:result.fraction_I_sensitivity},null,2));
