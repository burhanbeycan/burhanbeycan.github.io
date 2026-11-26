// node scripts/viscosity-gp-sensitivity.cjs [web-directory] [json-output]
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const defaultWeb = fs.existsSync(path.join(root, 'dist', 'viscosity-gp.js')) ? path.join(root, 'dist') : root;
const web = process.argv[2] ? path.resolve(process.argv[2]) : defaultWeb;
const output = process.argv[3] ? path.resolve(process.argv[3]) : path.join(web, 'viscosity-inputs', 'viscosity-gp-sensitivity.json');
const context = {window: {}};
for (const file of ['viscosity-engine.js', 'viscosity-intrinsic-data.js', 'viscosity-finite-data.js', 'viscosity-gp.js']) {
  vm.runInNewContext(fs.readFileSync(path.join(web, file), 'utf8'), context, {filename: file});
}
const {ViscosityGP: GP, VIS_FINITE, VIS_INTRINSIC} = context.window;
const grid = {length_u: [0.25, 0.5, 1], amplitudeCoefficient: [0.125, 0.25, 0.5], regularizationSD_eta_rel: 0.02};
const cases = {};
for (const length_u of grid.length_u) {
  for (const amplitudeCoefficient of grid.amplitudeCoefficient) {
    const config = {length_u, amplitudeCoefficient, regularizationSD_eta_rel: grid.regularizationSD_eta_rel};
    const result = GP.buildBenchmark(VIS_FINITE.series, VIS_INTRINSIC.rows, config);
    const chosen = [...result.within.filter(x => x.strategy === 'highest'), ...result.transfer];
    for (const trial of chosen) {
      const id = trial.seriesId + '/' + trial.strategy;
      cases[id] ??= {seriesId: trial.seriesId, strategy: trial.strategy, n: trial.folds.length,
        trainIds: trial.trainIds, testIds: trial.testIds, baseline: trial.metrics.baseline, grid: []};
      cases[id].grid.push({length_u, amplitudeCoefficient, gp: trial.metrics.gp});
    }
  }
}
for (const item of Object.values(cases)) {
  const values = item.grid.map(r => r.gp.rmse_pct);
  item.range_rmse_pct = [Math.min(...values), Math.max(...values)];
  item.nSmallerError = values.filter(x => x < item.baseline.rmse_pct).length;
}
const report = {protocol: 'Post-hoc, exploratory sensitivity on the same frozen splits. All nine equally reported combinations were specified before executing this grid. Validation, training data and fixed baseline are identical across settings. No setting selection or update to the reference GP.',
  grid, reference: GP.settings, cases};
fs.mkdirSync(path.dirname(output), {recursive: true});
fs.writeFileSync(output, JSON.stringify(report, null, 2) + '\n');
for (const [id, item] of Object.entries(cases)) console.log(id, 'Huggins=' + item.baseline.rmse_pct.toFixed(4),
  'GP range=' + item.range_rmse_pct.map(x => x.toFixed(4)).join('..'), 'smaller=' + item.nSmallerError + '/9');
