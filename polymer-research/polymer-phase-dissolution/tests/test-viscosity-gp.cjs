// node scripts/test-viscosity-gp.cjs [web-directory]
'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const fallback = fs.existsSync(path.join(root, 'dist', 'viscosity-gp.js')) ? path.join(root, 'dist') : root;
const web = process.argv[2] ? path.resolve(process.argv[2]) : fallback;
const context = {window: {}};
for (const file of ['viscosity-engine.js', 'viscosity-intrinsic-data.js', 'viscosity-finite-data.js', 'viscosity-gp.js']) {
  vm.runInNewContext(fs.readFileSync(path.join(web, file), 'utf8'), context, {filename: file});
}
const {ViscosityGP: GP, Viscosity: core, VIS_FINITE, VIS_INTRINSIC} = context.window;
const series = VIS_FINITE.series;
const intrinsic = VIS_INTRINSIC.rows;
function near(x, y, tol = 1e-9) { assert.ok(Math.abs(x-y) < tol * Math.max(1,Math.abs(x),Math.abs(y)), x + ' versus ' + y); }
const report = GP.buildBenchmark(series,intrinsic);
assert.equal(report.within.length, 12);
assert.equal(report.transfer.length, 2);
for (const check of [...report.within, ...report.transfer]) {
  for (const key of check.testIds) assert.ok(!check.trainIds.includes(key));
  for (const row of check.folds) {
    assert.ok(row.gpConditionalSD >= 0 && Number.isFinite(row.gpConditionalSD));
    assert.ok(row.gpConditionalSD <= GP.settings.amplitudeCoefficient * row.u ** 2 + 1e-12);
    assert.ok(row.physicallyAdmissible);
  }
}
for (const specimen of series) {
  for (const strategy of Object.keys(GP.strategiesTR)) {
    const fitted = GP.fitWithin(specimen,strategy);
    const out = GP.predict(fitted,0);
    near(out.gp, 1); near(out.gpConditionalSD,0); near(out.slope_dEtaRel_du, 1);
    near((GP.predict(fitted,1e-5).gp - 1) / (1e-5 / 1000), fitted.eta_mL_g, 1e-5);
    const bench = report.within.find(b => b.seriesId === specimen.id && b.strategy === strategy);
    const conventional = core.finiteBenchmark(specimen.rows);
    const conventionalFolds = strategy === 'highest' ? conventional.folds : conventional[strategy].folds;
    const altered = {...specimen, rows: specimen.rows.map(r => bench.testIds.includes(specimen.id+'/'+r.id) ? {...r,eta_rel:2*r.eta_rel} : r)};
    const copy = GP.fitWithin(altered, strategy);
    for (const row of bench.folds) {
      const actual = GP.predict(copy,row.c_g_L);
      near(actual.gp, row.gp);
      near(actual.baseline, conventionalFolds.find(r => r.id === row.id).huggins_eta_rel);
      if (strategy === 'lowInterpolation') near(actual.gp,actual.baseline);
    }
    assert.throws(() => GP.predict(fitted,-1));
  }
}
const theta = series.filter(s => s.solvent === 'Siklohekzan');
for (const sample of theta) {
  const train = theta.find(v => v.id !== sample.id);
  const original = GP.fitTransfer(train,sample,intrinsic);
  assert.equal(original.anchorIds.length,6);
  assert.equal(original.etaPriorDomain,'extrapolation');
  const blind = {...sample, sourceIntrinsic_L_g:10000};
  // Training for unseen grades must not even read the held-out viscosity rows.
  Object.defineProperty(blind, 'rows', {get() { throw Error('held-out targets accessed'); }});
  const clean = GP.fitTransfer(train, blind, intrinsic);
  for (const r of sample.rows) near(GP.predict(original,r.c_g_L).gp,GP.predict(clean,r.c_g_L).gp);
  const otherSolvent = {...sample, solvent: 'Toluen'};
  assert.throws(() => GP.fitTransfer(train,otherSolvent,intrinsic));
}
// Numerically independent reference posterior, computed by numpy.linalg.solve in ηrel units.
const copenhaver = report.within.find(r => r.seriesId === 'copenhaver-ps-toluene-30' && r.strategy === 'highest');
near(copenhaver.metrics.baseline.rmse_pct,7.405553764737902);
near(copenhaver.folds[0].gp,3.258154090834799);
near(copenhaver.folds[0].gpConditionalSD,0.11334753210089002);
assert.ok(copenhaver.metrics.gp.rmse > copenhaver.metrics.baseline.rmse);
const referencePrediction = copenhaver.folds[0].gp;
const sampleJCE = series.find(s => s.id === 'copenhaver-ps-toluene-30');
const settings = {...GP.settings, length_u: 1, amplitudeCoefficient: .125};
const alteredFit = GP.fitWithin(sampleJCE, 'highest', settings);
const distinctPrediction = GP.predict(alteredFit,19.93).gp;
assert.ok(Math.abs(distinctPrediction - referencePrediction) > .001);
settings.length_u = 3;
near(GP.predict(alteredFit,19.93).gp, distinctPrediction);
near(GP.predict(GP.fitWithin(sampleJCE, 'highest'),19.93).gp, referencePrediction);
assert.throws(() => GP.fitWithin(sampleJCE,'highest',{...GP.settings,length_u:0}));
const known = report.transfer.find(r => r.seriesId.endsWith('1770000'));
near(known.metrics.gp.rmse_pct,0.382,0.001);
console.log('PASS 4 measured series, 12 within-series folds, 2 whole-grade folds, withheld target mutation, disjoint IDs, physical limits, concentration units, covariance bounds and exact information parity.');
