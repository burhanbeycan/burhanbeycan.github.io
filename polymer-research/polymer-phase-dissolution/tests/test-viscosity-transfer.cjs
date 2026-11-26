'use strict';
const assert=require('node:assert/strict'), fs=require('node:fs'), path=require('node:path'), vm=require('node:vm');
const root=path.resolve(__dirname,'..');
const dir=fs.existsSync(path.join(root,'dist','viscosity-transfer-engine.js'))?path.join(root,'dist'):root;
const ctx={window:{}};vm.runInNewContext(fs.readFileSync(path.join(dir,'viscosity-intrinsic-data.js'),'utf8'),ctx);
const {makePairs,benchmark}=require(path.join(dir,'viscosity-transfer-engine.js'));
const pairs=makePairs(ctx.window.VIS_INTRINSIC.rows);
assert.equal(pairs.length,6);
for (const split of ['fraction','parent']){
  const f=benchmark(pairs,'cyclohexane',split), b=benchmark(pairs,'benzene',split);
  for(const [i,r] of f.entries()){
    assert.ok(!r.train.includes(r.id));
    if(split==='parent') assert.ok(r.train.every(id=>!id.startsWith(r.parent)));
    assert.ok(Math.abs(r.factor*b[i].factor-1)<1e-12);
  }
  const top=f.find(r=>r.id==='IK1500-1');
  assert.equal(top.observed,5306);
  assert.ok(Math.abs(top.unscaled-667*(506/228)**3)<1e-9);
  const altered=pairs.map(p=>p.id===top.id?{...p,benzene:{...p.benzene,eta_mL_g:3000}}:p);
  assert.equal(benchmark(altered,'cyclohexane',split).find(r=>r.id===top.id).adjusted,top.adjusted);
  console.log(split,top, 'MAPE', ...['baselineError','adjustedError'].map(key=>f.reduce((n,r)=>n+Math.abs(r[key]),0)/f.length));
}
