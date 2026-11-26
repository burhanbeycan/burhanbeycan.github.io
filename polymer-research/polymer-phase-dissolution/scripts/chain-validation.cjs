// node scripts/chain-validation.cjs [web-directory] [output-json]
// With no output file, write the full result to stdout.
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root=path.resolve(__dirname,'..');
const web=path.resolve(process.argv[2] || (fs.existsSync(path.join(root,'dist'))?path.join(root,'dist'):root));
const engine=require(path.join(web,'chain-engine.js'));
const bundle={window:{}};
for (const filename of ['chain-data.js','chain-external.js']) vm.runInNewContext(fs.readFileSync(path.join(web,filename),'utf8'),bundle,{filename});
const calibration=bundle.window.CHAIN_DATA, test=bundle.window.CHAIN_EXTERNAL;
// Two original material lots in Miyaki et al.; this is a same-publication
// sensitivity check. Both folds also extrapolate beyond the remaining masses.
const parentGroups = {
  BK2500: ['BK2500-4','BK2500-3','BK2500-2','BK2500-1'],
  IK1500: ['IK1500-2','IK1500-1']
};
const withinStudyByParent=Object.fromEntries(Object.entries(parentGroups).map(([parent,ids])=>{
  const kept=calibration.rows.filter(r=>!ids.includes(r.id));
  const withheld=calibration.rows.filter(r=>ids.includes(r.id));
  if(withheld.length!==ids.length)throw new Error('Ana numune kimlikleri eksik.');
  return [parent,{fitOn:kept.map(r=>r.id),nu:engine.fit(kept).nu,...engine.external(kept,withheld)}];
}));
const output={
  calibrationDOI:calibration.source.doi,externalDOI:test.source.doi,
  condition:{polymer:calibration.polymer,solvent:calibration.solvent,T_K:calibration.T_K},
  freeFit:engine.fit(calibration.rows),thetaFit:engine.fit(calibration.rows,.5),
  withinStudy:engine.benchmark(calibration.rows),withinStudyByParent,externalStudy:engine.external(calibration.rows,test.rows),
  measured:calibration.rows.concat(test.rows).map(r=>({
    id:r.id,Mw_g_mol:r.Mw_g_mol,Rg_nm:r.Rg_nm, cstar_g_L:engine.overlap(r.Mw_g_mol,r.Rg_nm,0).cstar_g_L,
    radiusSensitivity_g_L:r.Rg_reported_pm_nm==null?null:engine.overlapRadiusSensitivity(r.Mw_g_mol,r.Rg_nm,r.Rg_reported_pm_nm)
  }))
};
const serialized=JSON.stringify(output,null,2)+'\n';
if(process.argv[3])fs.writeFileSync(process.argv[3],serialized);else process.stdout.write(serialized);
