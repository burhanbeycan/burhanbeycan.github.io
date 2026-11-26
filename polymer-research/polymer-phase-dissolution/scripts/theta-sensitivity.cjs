'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'..');
const web=process.argv[2]?path.resolve(process.argv[2]):fs.existsSync(path.join(root,'dist','family-engine.js'))?path.join(root,'dist'):root;
const destination=process.argv[3]?path.resolve(process.argv[3]):path.join(web,'family-inputs','theta-sensitivity.json');
const M=require(path.join(web,'family-engine.js'));
const ctx={window:{}};vm.runInNewContext(fs.readFileSync(path.join(web,'family-intrinsic-data.js'),'utf8'),ctx);
const rows=ctx.window.THETA_FAMILY_DATA.rows, base=M.buildBenchmark(rows);
function compact(d) {
 return {n:d.design.nRows,familyCounts:Object.fromEntries(d.families.map(f=>[f.family,f.folds.length])),familyMAPE:Object.fromEntries(d.families.map(f=>[f.family,[f.scores.theta.mape,f.scores.ridge.mape]])),macro:d.overall.macro,micro:d.overall.micro};
}
// Intersect each family's observed mass interval. Cut determined by inputs (mass), never by targets.
const ranges=base.families.map(f=>[Math.min(...f.folds.map(v=>v.mass)),Math.max(...f.folds.map(v=>v.mass))]);
const massMin=Math.max(...ranges.map(r=>r[0])),massMax=Math.min(...ranges.map(r=>r[1]));
if(!(massMax>massMin))throw new Error('Ortak kütle aralığı boş.');
const overlap=rows.filter(r=>r.mass_g_mol>=massMin&&r.mass_g_mol<=massMax);
const precise=rows.filter(r=>r.temperature_basis!=='approximate_compilation');
const onlyMw=rows.filter(r=>r.mass_basis==='Mw');
const groups=[...new Set(onlyMw.map(r=>r.family))];
const report={
 target:'İntrinsik viskozite; tüm uyum ve sınamalar seçili altkümede baştan hesaplanır.',
 status:'Ana sonuçlar görüldükten sonra tanımlanan keşifsel duyarlılık sınamaları. Varsayılan model veya ceza ayarı değişmedi.',
 full:compact(base),
 excludeApproximateThetaTemperature:compact(M.buildBenchmark(precise)),
 commonMassRange:{min_g_mol:massMin,max_g_mol:massMax,...compact(M.buildBenchmark(overlap))},
 onlyExplicitMw:{n:onlyMw.length,families:groups,eligibleForThreeFoldProtocol:groups.length>=3,reason:'Polidien kaynağı mutlak M bildirir. Sadece özgün M_w etiketleri tutulursa üç aileli sınama yapılamaz.'},
};
fs.mkdirSync(path.dirname(destination),{recursive:true});fs.writeFileSync(destination,JSON.stringify(report,null,2)+'\n');
console.log('full',report.full.n,report.full.macro.theta.mape,report.full.macro.ridge.mape,
 'exact listed temperature',precise.length,report.excludeApproximateThetaTemperature.macro.theta.mape,report.excludeApproximateThetaTemperature.macro.ridge.mape,
 'common range',massMin,massMax,overlap.length,report.commonMassRange.macro.theta.mape,report.commonMassRange.macro.ridge.mape);
