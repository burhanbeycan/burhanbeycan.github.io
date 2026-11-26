'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const parent=path.resolve(__dirname,'..');
function base(name){return [process.argv[2],path.join(parent,'dist'),parent,__dirname,process.argv[3]].find(p=>p&&fs.existsSync(path.join(p,name)));}
const oclPath='chem/ocl-vendor.js', dir=base('repeat-match.js'), lib=base(oclPath), datasetPath=process.argv[4]||base('family-intrinsic-data.js'),enginePath=process.argv[5]||base('family-engine.js');
require(path.resolve(lib,oclPath));
const ocl=globalThis.PolyOcl(), reader=require(path.resolve(dir,'repeat-match.js')),catalog=require(path.resolve(datasetPath,'family-intrinsic-data.js')),models=require(path.resolve(enginePath,'family-engine.js'));
const keys=reader.registry(catalog,ocl),match=s=>reader.match(s,ocl,keys);
for(const [s,id] of [
  ['[*:1]CC(c1ccccc1)[*:2]','PS'],['c1cccc(c1)C([*:8])C[*:3]','PS'],['*C(C*)C1=CC=CC=C1','PS'],
  ['*CC(C)(c1ccccc1)*','PAMS'],['c1ccccc1C(C)(*)C*','PAMS'],
  ['*CC(C)(C)*','PIB'],['CC(*)(C)C*','PIB'],
  ['*CC=CC*','PBD'],['*CCC=C*','PBD'],
  ['*CC(C)=CC*','PI'],['*CC=C(C)C*','PI'],['*CCC=C(C)*','PI'],
]){
  const r=match(s);assert.equal(r.chemistry,id,s);
  const svg=reader.depiction(r);assert.ok(svg.includes('<svg')&&!/>Cl</.test(svg));
  assert.equal((svg.match(/>\*<\/text>/g)||[]).length,2);
}
for(const s of ['*CC(O)*','*CCC*','*CC(Cl)*','*CC(Br)*','*CC(*)C(C)C','*CC(*)C=C','*C[C@H](*)c1ccccc1','*C/C=C\\C*',
  '*C/C=C/C*','*C=C=C*','*CC(c1ccccc1)','[*:1]CC(c1ccccc1)[*:1]','*CC(*)c1ccccc1(','*CC(*)c1ccccc1-',
  '*CC(*)c1ccccc1()','**','*=CC=*','*CC.*','*CC(*)c1ccccc1[]','*CC(*)[13C]','*C[NH3+]C*','~CC*','*CC*.*CC*']){
  assert.throws(()=>match(s),undefined,s);
}
const selected=catalog.rows.find(r=>r.chemistry==='PIB'),model=models.fitForFamily(catalog.rows,selected.family),out=models.predict(model,selected.mass_g_mol,'PIB');
assert.ok(model.trainingIds.length>0&&!model.trainingIds.includes(selected.id));
const benchmark=models.buildBenchmark(catalog.rows).families.find(f=>f.family===selected.family).folds.find(f=>f.id===selected.id);
assert.equal(out.theta,benchmark.theta);assert.equal(out.ridge,benchmark.ridge);
console.log('PASS five chemistries, both anchor orientations, arbitrary atom order, two diene cuts, explicit aromatic Kekulé, blocked stereo/unsupported/invalid connectivity, attachment-only depiction and held-family predictions.');
