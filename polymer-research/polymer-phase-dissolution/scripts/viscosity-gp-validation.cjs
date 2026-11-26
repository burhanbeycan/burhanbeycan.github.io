// node scripts/viscosity-gp-validation.cjs [web-directory] [json-output]
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const defaultWeb = fs.existsSync(path.join(root, 'dist', 'viscosity-gp.js')) ? path.join(root, 'dist') : root;
const web = process.argv[2] ? path.resolve(process.argv[2]) : defaultWeb;
const output = process.argv[3] ? path.resolve(process.argv[3]) : path.join(web, 'viscosity-inputs', 'viscosity-gp-validation.json');
const context = {window: {}};
for (const file of ['viscosity-engine.js', 'viscosity-intrinsic-data.js', 'viscosity-finite-data.js', 'viscosity-gp.js']) {
  vm.runInNewContext(fs.readFileSync(path.join(web, file), 'utf8'), context, {filename: file});
}
const {VIS_FINITE, VIS_INTRINSIC, ViscosityGP} = context.window;
const out = ViscosityGP.buildBenchmark(VIS_FINITE.series, VIS_INTRINSIC.rows);
fs.mkdirSync(path.dirname(output), {recursive: true});
fs.writeFileSync(output, JSON.stringify(out, null, 2) + '\n');
for (const item of [...out.within, ...out.transfer]) {
  console.log(item.seriesId, item.strategy, 'n=' + item.folds.length, 'baseline RMSE%=' + item.metrics.baseline.rmse_pct.toFixed(3),
    'GP RMSE%=' + item.metrics.gp.rmse_pct.toFixed(3));
}
