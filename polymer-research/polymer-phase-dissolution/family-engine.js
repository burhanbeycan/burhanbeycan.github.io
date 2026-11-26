// Transfer across backbone groups at reported theta conditions, using only training labels for fitting.
(function (root) {
  'use strict';
  const referenceMass = 100000; // g/mol
  const penalty = 1; // Training weights sum to 1. After scaling, the unpenalized slope is halved.
  const atomicMass = {C:12.011, H:1.008};
  const chemistries = {
    PBD: {name:'1,4-polibütadien', family:'polydiene', formula:'C4H6', C:4, H:6, skeletalBonds:4},
    PI: {name:'1,4-poliizopren', family:'polydiene', formula:'C5H8', C:5, H:8, skeletalBonds:4},
    PIB: {name:'Poliizobütilen', family:'polyolefin', formula:'C4H8', C:4, H:8, skeletalBonds:2},
    PS: {name:'Polistiren', family:'vinylaromatic', formula:'C8H8', C:8, H:8, skeletalBonds:2},
    PAMS: {name:'Poli(α-metilstiren)', family:'vinylaromatic', formula:'C9H10', C:9, H:10, skeletalBonds:2},
  };
  const labels = {polydiene:'Polidienler',polyolefin:'Alifatik poliolefin (PIB)',vinylaromatic:'Vinil aromatikler'};
  function descriptor(id) {
    const mol = chemistries[id];
    if (!mol) throw new RangeError('Bu kimya için tekrarlayan birim tanımı yok.');
    const M0 = mol.C*atomicMass.C + mol.H*atomicMass.H;
    // Ratio divided by the fixed reference 1 mol/g to obtain a dimensionless number.
    return {formula:mol.formula, skeletalBonds:mol.skeletalBonds, M0, logBondsPerMass:Math.log(mol.skeletalBonds/M0)};
  }
  function groupBy(rows, keyOf) {
    const out = new Map();
    for (const row of rows) {
      const key = keyOf(row);
      if (!out.has(key)) out.set(key,[]);
      out.get(key).push(row);
    }
    return out;
  }
  function check(rows) {
    const ids = new Set();
    for (const r of rows) {
      if (!r.id || ids.has(r.id)) throw new Error('Ölçüm kimliği tekil olmalı.');
      ids.add(r.id);
      if (!chemistries[r.chemistry] || chemistries[r.chemistry].family !== r.family)
        throw new RangeError('Kimya ile aile eşleşmiyor.');
      if (![r.mass_g_mol,r.eta_mL_g].every(v=>Number.isFinite(v)&&v>0)) throw new RangeError('Ölçülmüş pozitif M ve [η] gerekli.');
      if (r.regime && !['theta','θ'].includes(r.regime)) throw new RangeError('Bu model θ kohortuna özgüdür.');
    }
  }
  function trainingRecords(rows) {
    const records=[];
    const families=groupBy(rows,r=>r.family);
    // Give each backbone group, each chemistry within it and each cited source within a chemistry equal weight.
    for (const familyRows of families.values()) {
      const species=groupBy(familyRows,r=>r.chemistry);
      for (const sameSpecies of species.values()) {
        const sources=groupBy(sameSpecies,r=>r.source_ref || r.ref || r.source_doi || 'unreported');
        for (const source of sources.values()) {
          const w=1/(families.size*species.size*sources.size*source.length);
          for (const r of source) records.push({
            w, x:descriptor(r.chemistry).logBondsPerMass,
            y:Math.log(r.eta_mL_g) - 0.5*Math.log(r.mass_g_mol/referenceMass),
          });
        }
      }
    }
    return records;
  }
  function fitForFamily(rows,heldFamily) {
    if (!labels[heldFamily]) throw new RangeError('Bilinmeyen aile.');
    // Do not access target viscosities here, even to validate the held-out rows.
    const train=rows.filter(r=>r.family!==heldFamily);
    check(train);
    const trainingFamilies=[...new Set(train.map(r=>r.family))].sort();
    if (trainingFamilies.length < 2) throw new RangeError('En az iki ayrı eğitim ailesi gerekli.');
    const a=trainingRecords(train);
    const weighted=fn=>a.reduce((s,t)=>s+t.w*fn(t),0);
    const intercept=weighted(t=>t.y), xMean=weighted(t=>t.x);
    const sd=Math.sqrt(weighted(t=>(t.x-xMean)**2));
    const slope=sd>1e-12? weighted(t=>(t.x-xMean)*(t.y-intercept)/sd)/(1+penalty):0;
    return {heldFamily,trainingFamilies,trainingIds:train.map(r=>r.id),nTrain:train.length,
      trainingChemistries:[...new Set(train.map(r=>r.chemistry))].sort(),intercept,xMean,xSD:sd,
      slope,massExponent:0.5,penalty,trainMassMin:Math.min(...train.map(r=>r.mass_g_mol)),trainMassMax:Math.max(...train.map(r=>r.mass_g_mol)),
      trainXMin:Math.min(...a.map(r=>r.x)),trainXMax:Math.max(...a.map(r=>r.x)),referenceMass,
    };
  }
  function predict(model,mass,chemistry) {
    if (!Number.isFinite(mass) || mass<=0) throw new RangeError('Mol kütlesi pozitif ve sonlu olmalı.');
    const {logBondsPerMass:x}=descriptor(chemistry);
    const logBaseline=model.intercept+0.5*Math.log(mass/model.referenceMass);
    const delta=model.xSD>1e-12?model.slope*(x-model.xMean)/model.xSD:0;
    const theta=Math.exp(logBaseline),ridge=Math.exp(logBaseline+delta);
    if (![theta,ridge].every(n=>n>0&&Number.isFinite(n))) throw new RangeError('Sayısal aralık dışında.');
    return {theta,ridge,domain:mass<model.trainMassMin||mass>model.trainMassMax?'extrapolation':'interpolation',
      descriptorDomain:x<model.trainXMin-1e-12||x>model.trainXMax+1e-12?'extrapolation':'interpolation'};
  }
  function metrics(rows,key,balanced=false) {
    const families=groupBy(rows,r=>r.family);
    const avg=fn=>{
      if (!balanced) return rows.reduce((s,r)=>s+fn(r),0)/rows.length;
      return [...families.values()].reduce((s,g)=>s+g.reduce((v,r)=>v+fn(r),0)/g.length,0)/families.size;
    };
    const d=r=>r[key]-r.observed, rel=r=>100*d(r)/r.observed;
    return {n:rows.length,mae:avg(r=>Math.abs(d(r))),rmse:Math.sqrt(avg(r=>d(r)**2)),mape:avg(r=>Math.abs(rel(r))),rmse_pct:Math.sqrt(avg(r=>rel(r)**2)),bias:avg(d)};
  }
  function scores(rows,balanced=false) {return {theta:metrics(rows,'theta',balanced),ridge:metrics(rows,'ridge',balanced)};}
  function buildBenchmark(rows) {
    check(rows);
    const familySets=groupBy(rows,r=>r.family);
    if(familySets.size<3) throw new RangeError('Bu sınama için üç farklı aile gerekli.');
    const families=[];
    for (const [family,test] of familySets) {
      const fit=fitForFamily(rows,family);
      const folds=test.map(r=>({id:r.id,chemistry:r.chemistry,family,mass:r.mass_g_mol,mass_basis:r.mass_basis,
        observed:r.eta_mL_g,ref:r.ref,source_ref:r.source_ref,source_doi:r.source_doi,solvent:r.solvent,T_C:r.T_C,temperature_basis:r.temperature_basis,
        ...predict(fit,r.mass_g_mol,r.chemistry)}));
      const sharedSources=[...new Set(test.map(r=>r.source_ref).filter(Boolean))]
        .filter(s=>rows.some(r=>r.family!==family&&r.source_ref===s));
      families.push({...fit,family,label:labels[family],folds,scores:scores(folds),
        support:{testChemistries:[...new Set(test.map(r=>r.chemistry))],sharedSourceRefs:sharedSources,
          outsideMass:folds.filter(r=>r.domain==='extrapolation').length,outsideDescriptor:folds.filter(r=>r.descriptorDomain==='extrapolation').length}});
    }
    families.sort((a,b)=>a.family.localeCompare(b.family));
    const all=families.flatMap(f=>f.folds);
    return {design:{massExponent:0.5,referenceMass,penalty,nRows:all.length,nFamilies:families.length,
      trainingWeights:'Eşit aile, aile içinde eşit kimya, kimya içinde eşit kaynak ve kaynak içinde eşit gözlem.',
      metric:'MAPE = ortalama(100 × |tahmin − ölçüm| / ölçüm). Aile dengeli ölçütte her saklı aile eşit ağırlıktadır. RMSE ve MAE mL/g, bağıl RMSE ve MAPE yüzde cinsindendir.',
      scaling:'ln([η]/(1 mL/g)) = b + 0,5 ln(M/(100000 g/mol)) + β z. z, ln(n_omurga / (M₀/(1 g/mol))) için eğitimde standartlaştırılmış değerdir. Temel modelde β=0.',
      scope:'Her veri satırı, kimyaya özgü bir θ çözeltisine aittir. Tekrar birimi ideal ana izomeri temel alır. Zincir rijitliği, çözücü kimliği, sıcaklık ve mikroyapı etkileri doğrudan modellenmez.',
      descriptors:Object.fromEntries(Object.keys(chemistries).map(k=>[k,{...chemistries[k],...descriptor(k)}])),
    },families,overall:{macro:scores(all,true),micro:scores(all)}};
  }
  const api={descriptor,fitForFamily,predict,buildBenchmark,labels,chemistries};
  if (root)root.FamilyModels=api;
  if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof window==='undefined'?null:window);
