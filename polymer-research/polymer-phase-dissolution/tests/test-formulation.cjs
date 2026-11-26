const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const root = path.resolve(__dirname,"..");
const web = fs.existsSync(path.join(root,"dist","formulation-model.js"))?path.join(root,"dist"):root;
const model = require(path.join(web,"formulation-model.js"));
const engine = require(path.join(web,"formulation-engine.js"));
const raw = fs.readFileSync(path.join(web,"viscosity-inputs","hpam-akbari.csv"),"utf8").trim().split(/\r?\n/);
const runRows = raw.slice(1).map(line=>line.split(","));
assert.equal(runRows.length,20);
assert.equal(model.nGroups,9);
assert.equal(model.test.length,15);
assert.equal(model.nRuns,runRows.length);
assert.equal(model.test.reduce((n,x)=>n+x.n,0),20);
assert.equal(model.test.find(r=>r.n===6).x.join("/"),"2000/12.5/30");
assert.equal(new Set(model.test.flatMap(x=>x.ids)).size,20);
for (const r of model.test) {
  assert.ok(engine.supported(model,r.x),r.x.join("/"));
  const bounds=engine.shearInterval(model,r.x.slice(0,2));
  assert.ok(r.x[2]>=bounds[0]-1e-7 && r.x[2]<=bounds[1]+1e-7);
}
assert.ok(!engine.supported(model,[3300,22.25,56]));
assert.ok(!engine.supported(model,[Infinity,12.5,20]));
assert.equal(engine.predict(model,[3300,22.25,56],"ridge"),null);
assert.equal(engine.shearInterval(model,[3300,22.25]),null);
for (const k of ["ridge","baseline","quadratic"]) {
  const expected=Math.exp(model.coef[k][0]);
  assert.equal(engine.predict(model,model.reference,k),expected);
  const errors=model.test.map(r=>r[k]-r.y);
  const rmse=Math.sqrt(errors.reduce((s,e)=>s+e*e,0)/errors.length);
  assert.ok(Math.abs(rmse-model.scores[k].rmse)<1e-10);
}
assert.ok(model.folds.every(g=>g.quadraticRank===10 && g.testConditions+g.trainConditions===15));
assert.ok(model.scores.quadratic.rmse<model.scores.ridge.rmse && model.scores.ridge.rmse<model.scores.baseline.rmse);
assert.ok(Math.abs(model.scores.quadratic.rmse-0.684863)<1e-5);
console.log("PASS: 20 measured runs, 9 formulation folds, 15 distinct conditions; hull, inference, rank and residuals verified.");
