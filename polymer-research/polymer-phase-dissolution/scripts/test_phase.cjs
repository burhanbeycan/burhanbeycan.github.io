"use strict";
const assert = require("node:assert/strict");
const root = require("node:fs").existsSync(require("node:path").join(__dirname, "../dist/index.html")) ? "../dist/" : "../";
const fh = require(root + "phase-engine.js");
global.window = {};
require(root + "measurements.js");
const system=window.EVIDENCE.systems.ps;
const nearly=(observed,expected,tolerance)=>assert.ok(Math.abs(observed-expected)<=tolerance,`${observed} ≠ ${expected} within ${tolerance}`);
const e=fh.equilibrium(100,0.65), result=fh.classify(0.02,e);
assert.equal(result.count,2);assert.equal(result.label,"meta");
nearly(result.curvature,0.2204081632653,1e-12);
nearly(e.bin[0],.0097680307278,1e-12);nearly(e.bin[1],.2719683621054,1e-12);
nearly(e.spin[0],1/26,1e-13);nearly(e.spin[1],.2,1e-13);
nearly(result.richFraction,.03902347956,1e-11);
assert.equal(fh.classify(.12,e).label,"unstable");
assert.equal(fh.classify(e.spin[0],e).label,"spinodal");
assert.equal(fh.classify(e.bin[0],e).label,"boundary");
assert.equal(fh.classify(.008,e).label,"one");
assert.equal(fh.classify(.1,fh.equilibrium(100,.59)).label,"one");
assert.equal(fh.classify(e.pc,fh.equilibrium(100,.605)).label,"critical");
for(const [N,C,d] of [[100,.0,1e-6],[32000,.16,1e-4],[3000,0,.015],[33000,.16,.02]]){
 const c0=fh.criticalChi(N,C)+d,q=fh.equilibrium(N,c0,C,24),finer=fh.equilibrium(N,c0,C,192);
 const[a,b]=q.bin;assert.ok(a<q.spin[0] && q.spin[1]<b);
 nearly(a,finer.bin[0],2e-9);nearly(b,finer.bin[1],2e-9);
 nearly(fh.df(a,N,c0,C),fh.df(b,N,c0,C),2e-9);
 nearly(fh.solventPotential(a,N,c0,C),fh.solventPotential(b,N,c0,C),2e-10);
}
for (const x of [.02,.4,.9]) {
 const N=420,C=.16,c0=.8,h=1e-5;
 nearly((fh.df(x+h,N,c0,C)-fh.df(x-h,N,c0,C))/(2*h),fh.d2(x,N,c0,C),5e-7);
}
let checked=0,largestError=0;
for(const s of system.samples)for(const record of system.measurements.filter(r=>r.series===s.id))for(const m of system.fits){
 const T=record.predicted_T_K[m.id],q=fh.equilibrium(s.N,m.A+m.B/T,m.C);
 assert.ok(q.bin);
 const mismatch=Math.min(...q.bin.map(x=>Math.abs(x-record.phi_p)));
 largestError=Math.max(largestError,mismatch);assert.ok(mismatch<4e-5);checked++;
}
console.log(JSON.stringify({acceptance:"pass",fixed_chi_metastable_point:result,grid_refinement:"24 and 192",predicted_temperatures_checked:checked,max_abs_phi_disagreement:largestError},null,2));
