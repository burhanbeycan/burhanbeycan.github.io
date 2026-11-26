// Dilute-coil size, geometric overlap and illustrative ideal-chain configurations.
(function (root) {
  'use strict';
  const NA = 6.02214076e23; // mol^-1, exact SI definition
  function positive(x, label) {
    if (!Number.isFinite(x) || x <= 0) throw new RangeError(label + ' sonlu ve pozitif olmalı.');
    return x;
  }
  function validate(rows) {
    if (!Array.isArray(rows) || !rows.length) throw new RangeError('Numune listesi boş.');
    const ids = new Set();
    for (const r of rows) {
      if (typeof r.id !== 'string' || !r.id || ids.has(r.id)) throw new RangeError('Her fraksiyon benzersiz bir kimlik taşımalı.');
      ids.add(r.id);
      positive(r.Mw_g_mol, 'M𝓌');
      positive(r.Rg_nm, 'R𝓰');
    }
  }
  // Equal weight per fraction in ln(Rg/nm). Centering the mass coordinate
  // makes the reference length interpretable and improves numerical conditioning.
  function fit(rows, fixedNu = null) {
    validate(rows);
    if (fixedNu !== null && !Number.isFinite(fixedNu)) throw new RangeError('Üs sonlu olmalı.');
    let meanX=0, meanY=0;
    for (const r of rows) {
      meanX += Math.log(r.Mw_g_mol) / rows.length;
      meanY += Math.log(r.Rg_nm) / rows.length;
    }
    let xx=0, xy=0;
    for (const r of rows) {
      const dx=Math.log(r.Mw_g_mol)-meanX;
      xx+=dx*dx; xy+=dx*(Math.log(r.Rg_nm)-meanY);
    }
    const distinct=new Set(rows.map(r=>r.Mw_g_mol)).size;
    if (fixedNu===null && distinct<2) throw new RangeError('Serbest üs için en az iki ayrı M𝓌 gerekir.');
    return {nu:fixedNu===null ? xy/xx : fixedNu, Rg_ref_nm:Math.exp(meanY), Mref_g_mol:Math.exp(meanX),
      minMw:Math.min(...rows.map(r=>r.Mw_g_mol)), maxMw:Math.max(...rows.map(r=>r.Mw_g_mol)), n:rows.length};
  }
  function predict(Mw, fitted) {
    positive(Mw,'M𝓌');
    positive(fitted.Rg_ref_nm,'Referans R𝓰'); positive(fitted.Mref_g_mol,'Referans M𝓌');
    if (!Number.isFinite(fitted.nu)) throw new RangeError('Model üssü geçersiz.');
    return fitted.Rg_ref_nm*Math.exp(fitted.nu*(Math.log(Mw)-Math.log(fitted.Mref_g_mol)));
  }
  function score(folds, key) {
    let abs=0, sq=0, relativeSq=0, logSq=0;
    for(const r of folds) {
      const e=r[key]-r.observed_nm;
      abs+=Math.abs(e); sq+=e*e;
      relativeSq+=(e/r.observed_nm)**2;
      logSq+=Math.log(r[key]/r.observed_nm)**2;
    }
    const n=folds.length;
    return {n, mae_nm:abs/n,rmse_nm:Math.sqrt(sq/n),rmse_pct:100*Math.sqrt(relativeSq/n),rmse_log:Math.sqrt(logSq/n)};
  }
  function fold(r, training, free, theta) {
    return {id:r.id,Mw_g_mol:r.Mw_g_mol,observed_nm:r.Rg_nm,free_nm:predict(r.Mw_g_mol,free),
      theta_nm:predict(r.Mw_g_mol,theta),domain:r.Mw_g_mol<free.minMw || r.Mw_g_mol>free.maxMw ? 'extrapolation':'interpolation',
      train_ids:training.map(t=>t.id),training_min_Mw:free.minMw,training_max_Mw:free.maxMw};
  }
  function report(folds) {
    return {folds,metrics:{free:score(folds,'free_nm'),theta:score(folds,'theta_nm')}};
  }
  function benchmark(rows) {
    validate(rows);
    if (rows.length<3) throw new RangeError('Fraksiyon dışarıda bırakma sınaması en az üç fraksiyon gerektirir.');
    return report(rows.map(r=>{
      const training=rows.filter(t=>t.id!==r.id);
      return fold(r,training,fit(training),fit(training,.5));
    }));
  }
  function external(training,test) {
    validate(training); validate(test);
    const trained=new Set(training.map(r=>r.id));
    if (test.some(r=>trained.has(r.id))) throw new RangeError('Bir fraksiyon hem eğitimde hem dış sınamada yer alamaz.');
    const free=fit(training), theta=fit(training,.5);
    return report(test.map(r=>fold(r,training,free,theta)));
  }
  // The spherical-volume convention uses a dilute Rg and the supplied molar mass.
  // One nm^3 = 10^-24 litre. For a polydisperse specimen this operational
  // Mw/Rg,z^3 ratio cannot reconstruct its entire distribution of coil volumes.
  function overlap(Mw, Rg, c_g_L) {
    positive(Mw,'M𝓌'); positive(Rg,'R𝓰');
    if (!Number.isFinite(c_g_L) || c_g_L<0) throw new RangeError('Derişim negatif olamaz ve sonlu olmalı.');
    const cstar_g_L=3*Mw/(4*Math.PI*(NA*1e-24)*Rg**3);
    return {cstar_g_L,ratio:c_g_L/cstar_g_L};
  }
  // Endpoint sensitivity at fixed Mw, not a confidence interval or a propagated
  // statistical standard deviation. The source does not define its ± as 1σ.
  function overlapRadiusSensitivity(Mw, Rg, reported_pm_nm) {
    positive(Mw,'M𝓌'); positive(Rg,'R𝓰');
    if (!Number.isFinite(reported_pm_nm) || reported_pm_nm<0 || reported_pm_nm>=Rg)
      throw new RangeError('Yarıçap aralığı tanımsız.');
    return {low_g_L:overlap(Mw,Rg+reported_pm_nm,0).cstar_g_L,
      high_g_L:overlap(Mw,Rg-reported_pm_nm,0).cstar_g_L};
  }
  // Reproducible Box–Muller sampling with an open-interval uniform generator.
  function normalGenerator(seed) {
    let state=seed>>>0, spare=null;
    function uniform() {
      state=(state+0x6D2B79F5)>>>0;
      let z=Math.imul(state^(state>>>15),state|1);
      z^=z+Math.imul(z^(z>>>7),z|61);
      return (((z^(z>>>14))>>>0)+.5)/4294967296;
    }
    return () => {
      if(spare!==null) {const y=spare;spare=null;return y;}
      const r=Math.sqrt(-2*Math.log(uniform())), angle=2*Math.PI*uniform();
      spare=r*Math.sin(angle);
      return r*Math.cos(angle);
    };
  }
  // Ideal, isolated, equal-mass Gaussian bead walk. "steps" is numerical
  // resolution, not polymer DP. The first walk is rendered, the following
  // independent walks estimate ensemble RMS. No per-conformation rescaling.
  function gaussian(Rg_nm,seed=17,steps=64,replicas=384) {
    positive(Rg_nm,'Referans R𝓰');
    if (!Number.isSafeInteger(seed) || seed<0 || seed>0xffffffff ||
        !Number.isInteger(steps) || steps<1 || steps>8192 ||
        !Number.isInteger(replicas) || replicas<1 || replicas>20000 || steps*replicas>4000000) {
      throw new RangeError('Geçersiz örnekleme parametreleri.');
    }
    const n=steps, beads=n+1;
    // <Rg²> = <b²> n(n+2)/(6(n+1)) for n+1 equally weighted beads.
    const stepRms_nm=Rg_nm*Math.sqrt(6*beads/(n*(n+2)));
    const sigma=stepRms_nm/Math.sqrt(3), normal=normalGenerator(seed);
    let sumRg2=0,sumRee2=0,sampleRg_nm,sampleRee_nm;
    const points=[];
    for (let walk=0;walk<=replicas;walk++) {
      let x=0,y=0,z=0,mx=0,my=0,mz=0,gyrationSum=0;
      if(walk===0)points.push([0,0,0]);
      for(let k=1;k<=n;k++) {
        x+=sigma*normal();y+=sigma*normal();z+=sigma*normal();
        const dx=x-mx,dy=y-my,dz=z-mz;
        mx+=dx/(k+1);my+=dy/(k+1);mz+=dz/(k+1);
        gyrationSum+=dx*(x-mx)+dy*(y-my)+dz*(z-mz);
        if(walk===0)points.push([x,y,z]);
      }
      const rg2=gyrationSum/beads,ree2=x*x+y*y+z*z;
      if(walk===0) {
        sampleRg_nm=Math.sqrt(rg2);sampleRee_nm=Math.sqrt(ree2);
        for (const p of points){p[0]-=mx;p[1]-=my;p[2]-=mz;}
      } else {sumRg2+=rg2;sumRee2+=ree2;}
    }
    return {points,sampleRg_nm,sampleRee_nm,ensembleRmsRg_nm:Math.sqrt(sumRg2/replicas),
      ensembleRmsRee_nm:Math.sqrt(sumRee2/replicas),targetRg_nm:Rg_nm,
      targetRee_nm:Math.sqrt(n)*stepRms_nm,stepRms_nm,seed,steps,replicas};
  }
  const api={NA,fit,predict,benchmark,external,overlap,overlapRadiusSensitivity,gaussian};
  root.ChainPhysics=api;
  if (typeof module!=='undefined' && module.exports) module.exports=api;
})(typeof window==='undefined'?globalThis:window);
