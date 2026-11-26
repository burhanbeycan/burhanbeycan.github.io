"use strict";
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { spawnSync } = require("node:child_process");
const {predict} = require("./inverse-engine.js");

const state = JSON.parse(fs.readFileSync(path.join(__dirname, "model_state.json"), "utf8"));
const queries = [];
for (const wt of [3.9, 4, 5, 6, 7.5, 9, 10, 10.1]) {
  for (const kelvin of [293.15, 298, 298.15, 300, 308.15, 309]) {
    for (const rate of [1, 10, 30, 60, 95, 101]) queries.push([wt, kelvin, rate]);
  }
}
const python = [
  "import json,sys",
  "from pva_model import predict_one",
  "state=json.load(open('model_state.json'))",
  "out=[predict_one(state,*row) for row in json.load(sys.stdin)]",
  "print(json.dumps(out))"
].join("\n");
const run = spawnSync(process.env.PYTHON || "python", ["-c", python], {
  cwd: __dirname, input: JSON.stringify(queries), encoding: "utf8"
});
assert.equal(run.status, 0, run.stderr);
const expected = JSON.parse(run.stdout);
for (let i = 0; i < queries.length; i++) {
  const actual = predict(state, ...queries[i]);
  assert.equal(actual.reason, expected[i].domain);
  assert.equal(actual.supported, expected[i].eta_app_cp !== null);
  if (actual.supported) assert.ok(Math.abs(actual.eta_cP - expected[i].eta_app_cp) < 1e-8);
}
assert.equal(predict(state, 4, 298, 0).supported, false);
assert.equal(predict(state, NaN, 298, 10).supported, false);
console.log(`Verified ${queries.length} Python/JS predictions and applicability flags.`);
