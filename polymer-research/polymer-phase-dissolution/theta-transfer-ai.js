// Theta-condition transfer. Shared between the static app and the reproducible benchmark.
(function (root) {
  'use strict';
  const model = typeof module === 'object' && module.exports
    ? require('./family-engine.js') : root.FamilyModels;
  const referenceMass = 100000; // g/mol. Fixed before testing, used only to select a reference.
  const slopes = [0, 0.25, 0.5, 0.75, 1]; // Fractions of the unpenalized slope.
  function groups(rows, key) {
    const result = new Map();
    for (const row of rows) {
      const k = row[key];
      if (!result.has(k)) result.set(k, []);
      result.get(k).push(row);
    }
    return result;
  }
  function mean(values) {
    return values.reduce((sum, value) => sum + value, 0) / values.length;
  }
  function weights(rows) {
    const families = groups(rows, 'family'), out = [];
    for (const family of families.values()) {
      const chemistries = groups(family, 'chemistry');
      for (const chemistry of chemistries.values()) {
        const sources = groups(chemistry, 'source_ref');
        for (const study of sources.values()) {
          const weight = 1 / (families.size * chemistries.size * sources.size * study.length);
          for (const r of study) out.push({r, weight});
        }
      }
    }
    return out;
  }
  const x = r => model.descriptor(r.chemistry).logBondsPerMass;
  const y = r => Math.log(r.eta_mL_g) - 0.5 * Math.log(r.mass_g_mol/referenceMass);
  function fit(rows) {
    if (!rows.length) throw new Error('Eğitim verisi boş.');
    const wr = weights(rows);
    const sum = fn => wr.reduce((s, a) => s + a.weight * fn(a.r), 0);
    const mx = sum(x), my = sum(y);
    const variance = sum(r => (x(r) - mx) ** 2);
    return {b: my, center: mx, slope: variance > 1e-14 ? sum(r => (x(r) - mx) * (y(r) - my)) / variance : 0};
  }
  function logPrediction(fit, r, shrink) {
    return fit.b + 0.5 * Math.log(r.mass_g_mol/referenceMass) + shrink * fit.slope * (x(r) - fit.center);
  }
  function tune(train) {
    const chem = groups(train, 'chemistry');
    const inner = [...chem.entries()].map(([code, test]) => ({
      test,
      fit: fit(train.filter(r => r.chemistry !== code)),
    }));
    const candidates = slopes.map(shrink => {
      const errors = inner.flatMap(f => f.test.map(r => ({...r,
        loss: (Math.log(r.eta_mL_g) - logPrediction(f.fit, r, shrink)) ** 2})));
      const loss = weights(errors).reduce((s, a) => s + a.weight * a.r.loss, 0);
      return {shrink, loss};
    });
    // A tie retains stronger shrinkage, including when an inner training set has constant x.
    const selected = candidates.reduce((best, c) => c.loss < best.loss - 1e-12 ? c : best);
    return {selected: selected.shrink, candidates, chemistryCount: chem.size};
  }
  function agreement(train) {
    // Per-chemistry discrepancy from out-of-chemistry predictions of the OLD, frozen Ridge.
    const held = [...groups(train, 'chemistry').entries()];
    const records = [];
    for (const [code, points] of held) {
      const trained = fit(train.filter(r => r.chemistry !== code));
      for (const r of points) records.push({...r, residual: Math.log(r.eta_mL_g) - logPrediction(trained, r, 0.5)});
    }
    const w = weights(records);
    let meanSquareOffset = 0;
    let withinSeries = 0;
    for (const set of groups(w.map(a => ({...a.r,weight:a.weight})), 'chemistry').values()) {
      const sum = set.reduce((s, v) => s + v.weight, 0);
      const avg = set.reduce((s, v) => s + v.weight * v.residual, 0) / sum;
      meanSquareOffset += sum * avg ** 2;
      // Within-chemistry scatter includes study differences when a chemical has multiple studies.
      // Effective df corrects the weighted sample variance, conditional on this chemistry.
      const sumSquared = set.reduce((s, r) => s + (r.weight/sum) ** 2, 0);
      if (sumSquared < 1 - 1e-12) withinSeries += set.reduce((s, v) => s + v.weight * (v.residual - avg) ** 2, 0)/(1 - sumSquared);
    }
    // Moment-based credibility; between-chemistry squared bias and individual residual variance
    // are both estimated without any observation from the target family. This is a heuristic
    // point-estimate shrinkage rule, not a measurement-error model or a calibrated posterior.
    const retention = meanSquareOffset + withinSeries > 0 ? meanSquareOffset / (meanSquareOffset + withinSeries) : 0;
    return {meanSquareOffset, withinSeries, retention};
  }
  function makeModel(rows, family) {
    const train = rows.filter(r => r.family !== family);
    if (groups(train, 'family').size < 2 || groups(train, 'chemistry').size < 2) throw new Error('İki eğitim ailesi gerekiyor.');
    for (const r of train) {
      if (![r.mass_g_mol, r.eta_mL_g].every(v => Number.isFinite(v) && v>0)) throw new Error('Geçersiz eğitim ölçümü.');
    }
    const fixed = fit(train);
    return {family, trainedIds: train.map(r => r.id), fixed, tuning: tune(train), agreement: agreement(train)};
  }
  function compatible(a,b) {
    return a.chemistry === b.chemistry && a.family === b.family && a.solvent === b.solvent &&
      Number.isFinite(a.T_C) && a.T_C === b.T_C && a.mass_basis === b.mass_basis;
  }
  function chooseReference(target, pool) {
    // The identity of the reference depends only on metadata and mass, never its viscosity.
    const available = pool.filter(r => compatible(target, r) && Number.isFinite(r.mass_g_mol) && r.mass_g_mol>0);
    return available.sort((a,b) => Math.abs(Math.log(a.mass_g_mol/referenceMass))-Math.abs(Math.log(b.mass_g_mol/referenceMass)) || a.id.localeCompare(b.id))[0] || null;
  }
  function predict(trained, target, reference=null) {
    if (target.family !== trained.family || model.chemistries[target.chemistry]?.family !== trained.family) throw new Error('Saklı aile ve tekrar yapısı uyuşmuyor.');
    if (!(Number.isFinite(target.mass_g_mol) && target.mass_g_mol > 0)) throw new Error('Mol kütlesi geçersiz.');
    const logOld = logPrediction(trained.fixed,target,0.5);
    const result = {baseline: Math.exp(logOld), tuned: Math.exp(logPrediction(trained.fixed,target,trained.tuning.selected))};
    if (!reference) return result;
    if (reference.id === target.id && target.id) throw new Error('Referans fraksiyonu kendisiyle sınanamaz.');
    if (!compatible(target, reference) || ![reference.mass_g_mol,reference.eta_mL_g].every(v=>Number.isFinite(v)&&v>0))
      throw new Error('Referansın kimyası, çözücüsü, sıcaklığı ve kütle tanımı hedefle eşleşmelidir.');
    const discrepancy = Math.log(reference.eta_mL_g) - logPrediction(trained.fixed,reference,0.5);
    result.control = reference.eta_mL_g * Math.sqrt(target.mass_g_mol/reference.mass_g_mol);
    result.adapted = Math.exp(logOld + trained.agreement.retention * discrepancy);
    return result;
  }
  function predictFromLibrary(rows, target, referencePool = rows) {
    if (!rows.some(r => compatible(target, r))) {
      throw new Error('Girilen çözücü, sıcaklık veya mol kütlesi tanımı bu kimyanın kayıtlı θ koşuluyla eşleşmiyor.');
    }
    const trained = makeModel(rows, target.family);
    // For an identified reference specimen, its own observation cannot supply its prediction.
    const eligible = referencePool.filter(r => !target.id || r.id !== target.id);
    const reference = chooseReference(target, eligible);
    const predictions = predict(trained, target, reference);
    return {predictions, reference: reference ? {
      id: reference.id, sample:reference.sample, chemistry:reference.chemistry,
      solvent:reference.solvent, T_C:reference.T_C,mass_basis:reference.mass_basis,
      mass_g_mol:reference.mass_g_mol,eta_mL_g:reference.eta_mL_g,
      source_doi:reference.source_doi, source:reference.source
    } : null, retention:trained.agreement.retention};
  }
  function describe(r) {
    return {id:r.id, chemistry:r.chemistry,family:r.family, mass:r.mass_g_mol, observed:r.eta_mL_g,
      solvent:r.solvent,T_C:r.T_C,mass_basis:r.mass_basis, source_doi:r.source_doi};
  }
  function report(points, unit, trainDetails) {
    const keys = points.some(p => 'control' in p) ? ['baseline','control','candidate'] : ['baseline','candidate'];
    function metrics(set) {
      const result = {};
      for (const name of keys) result[name] = mean(set.map(p => 100 * Math.abs(p[name]/p.observed - 1)));
      return result;
    }
    const byFamily = [...groups(points,unit).entries()].map(([family, rows]) => ({family,n:rows.length,...metrics(rows)}));
    const scores = Object.fromEntries(keys.map(k => [k,mean(byFamily.map(s=>s[k]))]));
    return {n:points.length,averagingUnit:unit,scores,microMAPE:metrics(points),byFamily,points,trainDetails};
  }
  function benchmark(rows) {
    if (new Set(rows.map(r=>r.id)).size !== rows.length || rows.some(r=>!r.id)) throw new Error('Numune kimlikleri tekil olmalı.');
    const zeros=[], refs=[], sources=[], details={};
    for(const [family, candidates] of groups(rows,'family')) {
      const trained=makeModel(rows,family);
      details[family]=trained;
      for(const target of candidates) {
        const init=predict(trained,target);
        zeros.push({...describe(target),baseline:init.baseline,candidate:init.tuned});
        const ref=chooseReference(target,candidates);
        if(ref && target.id!==ref.id) {
          const out=predict(trained,target,ref);
          refs.push({...describe(target),baseline:out.baseline,control:out.control,candidate:out.adapted,
            referenceId:ref.id,referenceMass:ref.mass_g_mol,referenceEta:ref.eta_mL_g,referenceDOI:ref.source_doi,
            sameStudy:target.source_ref===ref.source_ref});
        }
        const otherSource=chooseReference(target,candidates.filter(r=>r.source_ref!==target.source_ref));
        if(otherSource) {
          const out=predict(trained,target,otherSource);
          sources.push({...describe(target),baseline:out.baseline,control:out.control,candidate:out.adapted,
            referenceId:otherSource.id,referenceMass:otherSource.mass_g_mol,referenceEta:otherSource.eta_mL_g,referenceDOI:otherSource.source_doi,
            source:target.source_ref});
        }
      }
    }
    return {zeroShot:report(zeros,'family',details),oneReference:report(refs,'family'),externalPS:report(sources,'source')};
  }
  const api={makeModel,chooseReference,predict,predictFromLibrary,benchmark};
  if (root) root.ThetaTransferAI=api;
  if (typeof module === 'object' && module.exports) module.exports=api;
})(typeof window==='undefined'?null:window);
