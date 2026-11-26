/* Conserved binary Flory–Huggins / Cahn–Hilliard dynamics on a periodic square.
 * No clipping or mass repairs. Each IMEX step is accepted only inside (0,1)
 * and after verifying the discrete total free energy. dx=M=1 in reduced units. */
(function (scope) {
  'use strict';
  const tau = 2 * Math.PI;
  const freeEnergy = (p, N, chi, chi1=0) => p * Math.log(p) / N + (1-p)*Math.log1p(-p) + (chi+chi1*p)*p*(1-p);
  const curvature = (p, N, chi, chi1=0) => 1/(N*p) + 1/(1-p) - 2*chi + chi1*(2-6*p);
  const derivative = (p, N, chi, chi1=0) => (Math.log(p)+1)/N - Math.log1p(-p) - 1 + chi*(1-2*p) + chi1*(2*p-3*p*p);

  function fftLine(real, imag, n, offset, stride, inverse) {
    for (let i=1,j=0; i<n; i++) {
      let bit=n>>1;
      while (j & bit) {j ^= bit; bit >>= 1;}
      j ^= bit;
      if (i>=j) continue;
      const a=offset+i*stride, b=offset+j*stride;
      const r=real[a], v=imag[a]; real[a]=real[b]; imag[a]=imag[b]; real[b]=r; imag[b]=v;
    }
    for (let len=2; len<=n; len*=2) {
      const theta=(inverse ? tau : -tau)/len;
      const rc=Math.cos(theta), ic=Math.sin(theta), half=len/2;
      for (let i=0; i<n; i+=len) {
        let wr=1, wi=0;
        for (let j=0; j<half; j++) {
          const a=offset+(i+j)*stride, b=a+half*stride;
          const tr=wr*real[b]-wi*imag[b], ti=wi*real[b]+wr*imag[b];
          real[b]=real[a]-tr; imag[b]=imag[a]-ti; real[a]+=tr; imag[a]+=ti;
          const next=wr*rc-wi*ic; wi=wi*rc+wr*ic; wr=next;
        }
      }
    }
    if (inverse) for (let i=0; i<n; i++) {real[offset+i*stride]/=n; imag[offset+i*stride]/=n;}
  }
  function fft2(real,imag,n,inverse=false) {
    for (let i=0;i<n;i++) fftLine(real,imag,n,i*n,1,inverse);
    for (let i=0;i<n;i++) fftLine(real,imag,n,i,n,inverse);
  }
  function random(seed) {
    let x=seed>>>0;
    return () => {x=(Math.imul(x,1664525)+1013904223)>>>0; return x/4294967296;}
  }
  function mean(a) {let s=0; for (const x of a) s+=x; return s/a.length;}
  function shellPowers(sim) {
    const shells=new Float64Array(1+sim.size*sim.size/2);
    for (let i=1;i<sim.real.length;i++) shells[sim.shell[i]]+=sim.real[i]**2+sim.imag[i]**2;
    return shells;
  }
  function setup({N=100,chi=.85,chi1=0,phi=.15,mode='quench',size=64,kappa=1,seed=26011,noise=.001,dt=.5}={}) {
    if (!Number.isFinite(N) || N<1 || N>1e7 || !Number.isFinite(chi) || chi< -3 || chi>5 ||
        !Number.isFinite(chi1) || Math.abs(chi1)>5 || !Number.isFinite(phi) || phi<.0001 || phi>.9999 || !['quench','slab'].includes(mode) ||
        ![16,32,64,128].includes(size) || !Number.isFinite(kappa) || kappa<=0 ||
        !Number.isFinite(dt) || dt<=0 || !Number.isFinite(noise) || noise<0)
      throw Error('Giriş fiziksel veya sayısal çalışma aralığının dışında.');
    const count=size*size, vector=()=>new Float64Array(count);
    const sim={N,chi,chi1,phi,mode,size,kappa,seed,dt,time:0,steps:0,field:vector(), real:vector(),imag:vector(),
      mr:vector(),mi:vector(),trialReal:vector(),trialImag:vector(),trialField:vector(),tmp:vector(),
      shell:new Uint16Array(count),q2:vector(),early:null,curvature:curvature(phi,N,chi,chi1)};
    const wave=tau/size, rand=random(seed);
    const width=size/4, center=(size-1)/2, delta=Math.min(.92*(1-phi),2.7*phi);
    for (let y=0;y<size;y++) for (let x=0;x<size;x++) {
      const i=x+y*size, kx=Math.min(x,size-x), ky=Math.min(y,size-y);
      sim.shell[i]=kx*kx+ky*ky; sim.q2[i]=sim.shell[i]*wave*wave;
      if (mode==='quench') sim.field[i]=Math.min(noise,.12*Math.min(phi,1-phi))*(2*rand()-1);
      else {
        const band=(Math.tanh((x-center+width/2)/1.15)-Math.tanh((x-center-width/2)/1.15))/2;
        sim.field[i]=delta*band/(1-width/size);
      }
    }
    const rawMean=mean(sim.field);
    for (let i=0;i<count;i++) sim.field[i]+=phi-rawMean;
    sim.real.set(sim.field); fft2(sim.real,sim.imag,size);
    sim.energy=energy(sim,sim.field,sim.real,sim.imag);
    sim.initialEnergy=sim.energy;
    sim.initialMean=mean(sim.field);
    if (mode==='slab') {
      const p=profile(sim), centerLeft=size/2-1;
      sim.solventIn=1-p[centerLeft]; sim.solventOut=1-p[0];
      sim.initialSlab=slabMetrics(sim,p);
    } else {
      sim.initialPowers=shellPowers(sim);
      sim.currentPrediction=dominantMode(sim);
    }
    return sim;
  }
  function energy(sim, field, real, imag) {
    let bulk=0, grad=0;
    for (let i=0;i<field.length;i++) {
      const v=field[i]; if (!(v>0 && v<1)) return Infinity;
      bulk+=freeEnergy(v,sim.N,sim.chi,sim.chi1);
      if (i) grad+=sim.q2[i]*(real[i]*real[i]+imag[i]*imag[i]);
    }
    return bulk/field.length + .5*sim.kappa*grad/(field.length*field.length);
  }
  function step(sim, maxDt=sim.dt) {
    const n=sim.field.length;
    if (!(maxDt>0 && Number.isFinite(maxDt))) throw Error('Geçersiz zaman adımı.');
    for (let i=0;i<n;i++) {sim.mr[i]=derivative(sim.field[i],sim.N,sim.chi,sim.chi1); sim.mi[i]=0;}
    fft2(sim.mr,sim.mi,sim.size);
    let dt=sim.lastDt ? Math.min(maxDt,1.4*sim.lastDt) : maxDt;
    for (let attempt=0;attempt<40;attempt++,dt*=.5) {
      for (let i=0;i<n;i++) {
        const q=sim.q2[i], scale=1/(1+dt*sim.kappa*q*q);
        sim.trialReal[i]=(sim.real[i]-dt*q*sim.mr[i])*scale;
        sim.trialImag[i]=(sim.imag[i]-dt*q*sim.mi[i])*scale;
      }
      sim.trialField.set(sim.trialReal); sim.tmp.set(sim.trialImag);
      fft2(sim.trialField,sim.tmp,sim.size,true);
      const candidate=energy(sim,sim.trialField,sim.trialReal,sim.trialImag);
      if (candidate>sim.energy+5e-14*(1+Math.abs(sim.energy))) continue;
      [sim.real,sim.trialReal]=[sim.trialReal,sim.real]; [sim.imag,sim.trialImag]=[sim.trialImag,sim.imag];
      [sim.field,sim.trialField]=[sim.trialField,sim.field];
      sim.energy=candidate; sim.time+=dt; sim.steps++; sim.lastDt=dt;
      if (sim.mode==='quench' && sim.currentPrediction && (!sim.early || sim.time<=2.5+1e-9))
        sim.early=measuredMode(sim);
      return dt;
    }
    throw Error('Kararlı zaman adımı bulunamadı. Son kabul edilen alan korunuyor.');
  }
  function dominantMode(sim) {
    const fpp=sim.curvature, wave=tau/sim.size;
    if (fpp>=0) return null;
    let best=0, peak=null;
    for (const k2 of sim.shell) {
      const q2=wave*wave*k2, sigma=-q2*(fpp+sim.kappa*q2);
      if (sigma>best) {best=sigma; peak={k2,q:Math.sqrt(q2),sigma,lambda:tau/Math.sqrt(q2)};}
    }
    return peak;
  }
  function measuredMode(sim) {
    const P=shellPowers(sim), wave=tau/sim.size, cutoff=-sim.curvature/sim.kappa;
    let best=-Infinity, peak=null, maxP=0;
    for(let j=1;j<P.length;j++) maxP=Math.max(maxP,sim.initialPowers[j]);
    for(let j=1;j<P.length;j++) {
      const q2=j*wave*wave;
      if (!(q2<cutoff) || sim.initialPowers[j]<=Math.max(1e-28,maxP*1e-8) || P[j]<=0) continue;
      const rate=Math.log(P[j]/sim.initialPowers[j])/(2*sim.time);
      if(rate>best) {best=rate; peak={k2:j,q:Math.sqrt(q2),lambda:tau/Math.sqrt(q2),sigma:rate,t:sim.time};}
    }
    return peak;
  }
  function profile(sim) {
    const n=sim.size, p=new Float64Array(n);
    for(let i=0;i<sim.field.length;i++) p[i%n]+=sim.field[i]/n;
    return p;
  }
  function frontAt(p, threshold) {
    let found=null;
    for(let x=0;x<p.length/2-1;x++) {
      const a=1-p[x], b=1-p[x+1];
      if ((a>=threshold && b<threshold) || (a<threshold && b>=threshold)) {
        if(found!==null) return null;
        found=x+(threshold-a)/(b-a);
      }
    }
    return found;
  }
  function slabMetrics(sim,p=profile(sim)) {
    const n=sim.size;
    let sum=0, count=0;
    for(let x=0;x<n;x++) if(Math.abs(x-(n-1)/2)<n/8) {sum+=1-p[x];count++;}
    const contrast=sim.solventOut-sim.solventIn;
    const cut=frac=>frontAt(p,sim.solventIn+frac*contrast);
    const x10=cut(.1), x90=cut(.9), x20=cut(.2);
    return {solvent:sum/count,x20,width:x10!==null&&x90!==null?Math.abs(x10-x90):null};
  }
  const api={setup,step,profile,slabMetrics,freeEnergy,derivative,curvature,dominantMode,mean,fft2};
  if (typeof module!=='undefined' && module.exports) module.exports=api;
  scope.PolyCH=api;
})(globalThis);
