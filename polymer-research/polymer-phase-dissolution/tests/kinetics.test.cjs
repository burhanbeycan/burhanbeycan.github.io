'use strict';
const assert=require('node:assert/strict');
const root = require('node:fs').existsSync(require('node:path').join(__dirname, '../dist/index.html')) ? '../dist/' : '../';
const CH = require(root + 'kinetics.js');

function near(a,b,tol,label){assert.ok(Math.abs(a-b)<tol,label+' '+a+' compared with '+b);}
function stableSteps(s,n,maxDt=s.dt){
  for(let i=0;i<n;i++){
    const e=s.energy,t=s.time;CH.step(s,maxDt);
    near(CH.mean(s.field),s.initialMean,2e-13,'mass');
    assert.ok(s.energy<=e+1e-13,'non-increasing free energy');
    assert.ok(s.time>t && s.lastDt<=maxDt+1e-15,'adaptive step');
    assert.ok(s.field.every(v=>v>0 && v<1),'logarithmic domain');
  }
}

const n=16;
const a=Float64Array.from({length:n*n},(_,i)=>Math.sin(i*.71)+.2*Math.cos(i*.19)), b=new Float64Array(n*n);
const a0=new Float64Array(a);
CH.fft2(a,b,n);CH.fft2(a,b,n,true);
near(Math.max(...a.map((v,i)=>Math.abs(v-a0[i]))),0,2e-13,'FFT round trip');
near(Math.max(...b.map(Math.abs)),0,2e-13,'FFT imaginary leakage');

for(const p of [.035,.31,.82]) for(const C of [-.25,0,.4]) {
  const delta=1e-5, c=.61,N=97.4;
  near((CH.freeEnergy(p+delta,N,c,C)-CH.freeEnergy(p-delta,N,c,C))/(2*delta),CH.derivative(p,N,c,C),1e-8,'potential derivative at '+p);
  near((CH.derivative(p+delta,N,c,C)-CH.derivative(p-delta,N,c,C))/(2*delta),CH.curvature(p,N,c,C),2e-8,'curvature at '+p);
}
near(CH.curvature(.02,100,.65),.220408163265306,1e-13,'metastable reference');

const u=CH.setup({noise:1e-6,dt:.05});
stableSteps(u,50);
assert.equal(u.currentPrediction.k2,25);
assert.equal(u.early.k2,25);
near(u.early.lambda,12.8,1e-10,'measured wavelength');
const q2=u.currentPrediction.q**2;
const sigmaStep=Math.log((1-u.dt*q2*u.curvature)/(1+u.kappa*u.dt*q2*q2))/u.dt;
near(u.early.sigma,sigmaStep,1e-5,'IMEX dispersion');
near(u.early.sigma,u.currentPrediction.sigma,.00035,'continuum dispersion');

const stable=CH.setup({phi:.02,N:100,chi:.65});
const rms=s=>Math.sqrt(s.field.reduce((a,v)=>a+(v-s.phi)**2,0)/s.field.length);
const rms0=rms(stable);
stableSteps(stable,25);assert.equal(stable.currentPrediction,null);assert.ok(rms(stable)<rms0);

const fine=CH.setup({size:64,N:1600,phi:.11,chi:.71,chi1:.13,noise:.0008});
stableSteps(fine,70);assert.ok(fine.energy<fine.initialEnergy);

const sheet=CH.setup({mode:'slab',N:100,phi:.2,chi:.30});
stableSteps(sheet,10);
let m=CH.slabMetrics(sheet);
assert.ok(m.x20>sheet.initialSlab.x20 && m.solvent>sheet.initialSlab.solvent && m.width>sheet.initialSlab.width,'penetration and swelling in fixed slab');
stableSteps(sheet,350);m=CH.slabMetrics(sheet);assert.equal(m.x20,null);assert.equal(m.width,null);
const stiff=CH.setup({mode:'slab',size:32,N:30000,phi:.9999,chi:.59,chi1:.17});
stableSteps(stiff,8);assert.ok(stiff.lastDt<.1,'backtracking near endpoint');
assert.throws(()=>CH.setup({phi:0}));assert.throws(()=>CH.setup({N:NaN}));assert.throws(()=>CH.setup({chi1:Infinity}));
console.log('PASS FFT; FH reduced free energy / chemical potential / curvature for three C values; mass and energy at EVERY step; FH reference curvature; positivity; two initial conditions; front detection; linear dispersion and wavelength; adaptive step at endpoint.');
console.log('Default dilute quench (noise=1e-6, dt*=0.05, t*=2.5): wavelength measured / discrete prediction',u.early.lambda,u.currentPrediction.lambda,'and growth rate measured / continuous-time discrete-space',u.early.sigma,u.currentPrediction.sigma);
