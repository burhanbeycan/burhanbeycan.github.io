'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const project = path.resolve(__dirname, '..');
const dist = path.resolve(process.argv[2] || (fs.existsSync(path.join(project, 'dist/activity-engine.js')) ? path.join(project, 'dist') : project));
global.window = global;
for (const file of ['phase-engine.js','osmotic-engine.js','dense-engine.js','dense-results.js','activity-engine.js']) {
  require(path.join(dist, file));
}
const activity = global.SolventActivity;
const close = (value, expected, tol = 1e-10) => {
  assert.ok(Math.abs(value - expected) < tol, value + ' vs ' + expected);
};
const csv = fs.readFileSync(path.join(dist, 'transfer-inputs/ps_cyclohexane_activity.csv'), 'utf8').trim().split(/\r?\n/);
assert.equal(csv.length, 6);
const rows = activity.compare();
assert.equal(rows.length, 5);
for (const [i, row] of rows.entries()) {
  const fields = csv[i + 1].split(',');
  assert.equal(fields[1], 'III');
  assert.equal(fields[3], 'II');
  assert.equal(fields[8], row.phiPrinted);
  assert.equal(fields[10], row.minusLogAPrinted);
  assert(row.phi > 259 * 0.934 / 1000);
  for (const model of ['observed','baseline','extended','guard']) {
    const r = row[model];
    close(-Math.log(r.activity), r.minusLogA);
    close(r.deltaMu_J_mol, 8.31446261815324 * row.T_K * r.lnA);
    if (model !== 'observed') {
      close(-r.deltaMu_J_mol, r.Pi_kPa * activity.referenceVolumeMlMol / 1000);
      assert.equal(r.phases, 1);
    }
  }
  close(row.extended.minusLogA - row.baseline.minusLogA,
    activity.D() * row.phi ** 3 * (2 - 3 * row.phi));
}
for (const model of ['baseline','extended','guard']) {
  const pure = activity.predict(307.15, 0, model);
  assert.equal(pure.activity, 1);
  const pinned = activity.predict(317.15, 2/3, model);
  close(pinned.activity, activity.predict(317.15, 2/3, 'baseline').activity);
}
const s = global.DENSE_RESULTS.sample;
const frozen = global.DenseFH.state(297.15, s.Mn_g_mol, s.vbar_ml_g, s.Vs_ml_mol, global.DENSE_RESULTS.base, activity.D());
assert(frozen.eq.bin);
const [lo, hi] = frozen.eq.bin;
const edge1 = activity.predict(297.15, lo + (hi-lo)/3, 'extended');
const edge2 = activity.predict(297.15, lo + 2*(hi-lo)/3, 'extended');
assert.equal(edge1.phases, 2);
assert.equal(edge2.phases, 2);
close(edge1.activity, edge2.activity);
const scores = activity.metrics();
close(scores.baseline.maeMinusLogA, 0.0174930745582, 1e-9);
close(scores.extended.maeMinusLogA, 0.0143507761737, 1e-8);
close(scores.guard.maeMinusLogA, 0.0162079473879, 1e-8);
assert.equal(activity.metrics([]).extended.maeMinusLogA, null);
assert.throws(() => activity.predict(0, 0.4), RangeError);
assert.throws(() => activity.predict(307.15, 1), RangeError);
assert.throws(() => activity.predict(307.15, 0.4, 'random'), RangeError);
console.log('PASS Table III experimental precision, activity/chemical-potential identities, pure-solvent reference, coexistence plateau, correction at 2/3, five held-out rows.');
