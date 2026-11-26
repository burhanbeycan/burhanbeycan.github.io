// node scripts/viscosity-validation.cjs [web-directory] [report.json]
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const defaultWeb = fs.existsSync(path.join(root, 'dist', 'viscosity-engine.js')) ? path.join(root, 'dist') : root;
const web = process.argv[2] ? path.resolve(process.argv[2]) : defaultWeb;
const output = process.argv[3] ? path.resolve(process.argv[3]) : path.join(web, 'viscosity-inputs', 'viscosity-validation.json');
const model = require(path.join(web, 'viscosity-engine.js'));

function loadData(file, key) {
  const context = {window: {}};
  vm.runInNewContext(fs.readFileSync(path.join(web, file), 'utf8'), context, {filename: file});
  return context.window[key];
}

const intrinsic = loadData('viscosity-intrinsic-data.js', 'VIS_INTRINSIC');
const external = loadData('viscosity-external-data.js', 'VIS_EXTERNAL');
const finite = loadData('viscosity-finite-data.js', 'VIS_FINITE');
const systems = {};
for (const key of new Set(intrinsic.rows.filter(r => r.training === true).map(r => r.system))) {
  const train = intrinsic.rows.filter(r => r.training === true && r.system === key);
  const extra = intrinsic.rows.filter(r => r.series === 'additional' && r.system === key);
  const other = external.rows.filter(r => r.system === key);
  const checks = model.studyBenchmark(train, other);
  const blocked = {};
  for (const parent of new Set(train.map(r => r.parent).filter(Boolean))) {
    const test = train.filter(r => r.parent === parent);
    const remainder = train.filter(r => r.parent !== parent);
    if (remainder.length >= 2) blocked[parent] = model.mhReport(model.mhSet(remainder, test));
  }
  systems[key] = {
    fit: model.fitMH(train), training_ids: train.map(r => r.id), loo: checks.loo,
    parentBlocked: blocked, extraSameStudy: extra.length ? model.mhReport(model.mhSet(train, extra)) : null,
    external: checks.external,
  };
}
const series = {};
for (const sample of finite.series) {
  const checks = model.finiteBenchmark(sample.rows);
  series[sample.id] = {
    ...checks, source: sample.source.url, temperature_K: sample.T_K, molarMass_g_mol: sample.Mw_g_mol,
    // The full-series intercept and the authors' value reuse observations. Neither enters the held-out scoring.
    fullSeriesFit: model.fitDilute(sample.rows),
    sourceIntrinsic_mL_g: sample.sourceIntrinsic_L_g == null ? null : sample.sourceIntrinsic_L_g * 1000,
  };
}
const report = {
  version: 1,
  units: {mass: 'g/mol', intrinsic: 'mL/g', concentration: 'g/L', relativeViscosity: 'dimensionless'},
  intrinsicProtocol: 'For each solvent, fit six paired fractions with equal-weight OLS in log coordinates. Use measured Mw as the empirical mass coordinate. Score other-study measurements at the same solvent and temperature, without using them during fitting. Parent groups and additional fractions from the training study are reported separately.',
  finiteProtocol: 'Per specimen, fit reduced and inherent viscosity lines on the lower concentrations. Hold out the largest distinct c, or two largest as a sensitivity check. For low-concentration interpolation, hold the middle of the three lowest levels and fit only its two neighbours. Retain flagged points. All series use the relative viscosities reported by the source.',
  scope: 'Huggins/Kraemer are dilute expansions. Use the reported u and residuals to interpret extrapolations. Two-point fits have zero residual degrees of freedom. Each mPa·s prediction requires the matched solvent dynamic viscosity.',
  systems, series,
};
fs.mkdirSync(path.dirname(output), {recursive: true});
fs.writeFileSync(output, JSON.stringify(report, null, 2) + '\n');
for (const [key, result] of Object.entries(systems)) {
  console.log(key, 'a=' + result.fit.a.toFixed(5), 'LOO RMSE%=' + result.loo.metrics.rmse_pct.toFixed(3),
    'external n=' + (result.external?.metrics.n ?? 0), 'external RMSE%=' + (result.external?.metrics.rmse_pct.toFixed(3) ?? '—'));
}
for (const [key, result] of Object.entries(series)) {
  console.log(key, 'highest-c Huggins RMSE%=' + result.metrics.huggins.rmse_pct.toFixed(3),
    'lowest-three interpolation RMSE%=' + result.lowInterpolation.metrics.huggins.rmse_pct.toFixed(3));
}
