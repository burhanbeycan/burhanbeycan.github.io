/* Huggins mean + Gaussian-process residual for measured dilute-solution series.
   Prediction uses only the training observations. All concentration inputs are g/L. */
(function (root) {
  'use strict';
  const core = root.Viscosity || (typeof require === 'function' ? require('./viscosity-engine.js') : null);
  if (!core) throw new Error('Önce viscosity-engine.js yüklenmeli.');
  // Fixed before evaluation; numerical regularization in ηr units is separate from reported experimental SDs.
  const settings = Object.freeze({length_u: 0.5, amplitudeCoefficient: 0.25, regularizationSD_eta_rel: 0.02});
  const strategies = ['highest', 'upperTwo', 'lowInterpolation'];
  const strategiesTR = {highest: 'En yüksek derişim', upperTwo: 'En yüksek iki derişim', lowInterpolation: 'Alt aralıkta ara değer'};
  const matchedSystem = 'ps-cyclohexane-theta';

  function validateConfig(config) {
    const out = {};
    for (const key of Object.keys(settings)) {
      if (!(Number.isFinite(config[key]) && config[key] > 0)) throw new RangeError('GP parametresi pozitif ve sonlu olmalı: ' + key);
      out[key] = config[key];
    }
    return Object.freeze(out);
  }
  function kernel(u, v, config) {
    const d = (u - v) / config.length_u;
    return config.amplitudeCoefficient ** 2 * u ** 2 * v ** 2 * Math.exp(-0.5 * d ** 2);
  }
  function dot(a, b) {
    return a.reduce((sum, value, i) => sum + value * b[i], 0);
  }
  // Cholesky solves for short positive-definite covariance matrices in this application.
  function factor(matrix) {
    const rows = matrix.map(row => new Array(row.length).fill(0));
    for (let i = 0; i < rows.length; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = matrix[i][j];
        for (let k = 0; k < j; k++) sum -= rows[i][k] * rows[j][k];
        if (i === j && !(sum > 0)) throw new RangeError('Kovaryans matrisi pozitif tanımlı olmalı.');
        rows[i][j] = i === j ? Math.sqrt(sum) : sum / rows[j][j];
      }
    }
    return rows;
  }
  function forward(lower, rhs) {
    const z = [];
    for (let i = 0; i < rhs.length; i++) {
      let sum = rhs[i];
      for (let j = 0; j < i; j++) sum -= lower[i][j] * z[j];
      z.push(sum / lower[i][i]);
    }
    return z;
  }
  function backwardTranspose(lower, rhs) {
    const x = new Array(rhs.length);
    for (let i = rhs.length - 1; i >= 0; i--) {
      let sum = rhs[i];
      for (let j = i + 1; j < x.length; j++) sum -= lower[j][i] * x[j];
      x[i] = sum / lower[i][i];
    }
    return x;
  }
  function baselineAt(u, kH) { return 1 + u + kH * u ** 2; }
  function trainResidual(sample, rows, eta, kH, config) {
    if (!Number.isFinite(kH)) throw new RangeError('Huggins katsayısı sonlu olmalı.');
    const uTrain = rows.map(r => core.scale(eta, r.c_g_L).u);
    const deviations = rows.map((r, i) => r.eta_rel - baselineAt(uTrain[i], kH));
    if (!deviations.every(Number.isFinite)) throw new RangeError('Uyum için sonlu yanıtlar gerekli.');
    const covariance = uTrain.map((u, i) => uTrain.map((v, j) =>
      kernel(u, v, config) + (i === j ? config.regularizationSD_eta_rel ** 2 : 0)));
    const lower = factor(covariance);
    const weights = backwardTranspose(lower, forward(lower, deviations));
    return {eta_mL_g: eta, kH, settings: config, uTrain, lower, weights, minU: Math.min(...uTrain), maxU: Math.max(...uTrain),
      trainSeriesId: sample.id, trainIds: rows.map(r => sample.id + '/' + r.id),
      trainingObservations: rows.map((r, i) => ({id: r.id, u: uTrain[i], c_g_L: r.c_g_L, observed: r.eta_rel,
        baseline: baselineAt(uTrain[i], kH), residual: deviations[i]})),
    };
  }
  function splitRows(sample, strategy) {
    if (!strategies.includes(strategy)) throw new RangeError('Geçersiz sınama seçimi.');
    const sorted = [...sample.rows].sort((a, b) => a.c_g_L - b.c_g_L);
    const levels = [...new Set(sorted.map(r => r.c_g_L))];
    if (new Set(sorted.map(r => r.id)).size !== sorted.length) throw new RangeError('Ölçüm kimlikleri benzersiz olmalı.');
    if (!levels.every(c => Number.isFinite(c) && c > 0)) throw new RangeError('Sınama için pozitif derişimler gerekli.');
    const interior = strategy === 'lowInterpolation';
    const count = strategy === 'upperTwo' ? 2 : 1;
    if (levels.length < count + 2) throw new RangeError('Yeterli farklı derişim yok.');
    const testLevels = new Set(interior ? [levels[1]] : levels.slice(-count));
    const used = interior ? sorted.filter(r => r.c_g_L <= levels[2]) : sorted;
    return {train: used.filter(r => !testLevels.has(r.c_g_L)), test: used.filter(r => testLevels.has(r.c_g_L))};
  }
  function fitWithin(sample, strategy = 'highest', config = settings) {
    const {train} = splitRows(sample, strategy);
    const fit = core.fitDilute(train);
    return {...trainResidual(sample, train, fit.eta_mL_g, fit.kH, validateConfig(config)), testSeriesId: sample.id, strategy,
      etaPriorSource: 'Huggins · yalnız eğitim derişimleri', etaPriorDomain: 'within-study-fit'};
  }
  function eligibleGrade(sample) {
    // No surrogate structure-to-solubility or temperature mapping. The benchmark has two eligible grades.
    return sample.polymer === 'Polistiren' && (sample.solvent === 'Siklohekzan' || sample.solvent === 'cyclohexane') &&
      Math.abs(sample.T_K - 307.65) < 1e-8 && Number.isFinite(sample.Mw_g_mol) && sample.Mw_g_mol > 0;
  }
  function fitTransfer(trainSample, testSample, intrinsicRows, config = settings) {
    if (trainSample.id === testSample.id || trainSample.Mw_g_mol === testSample.Mw_g_mol ||
        !eligibleGrade(trainSample) || !eligibleGrade(testSample)) throw new RangeError('Aktarım farklı mol kütlelerine ve aynı PS–siklohekzan koşuluna özgü.');
    const anchors = intrinsicRows.filter(r => r.system === matchedSystem && r.training === true && r.study === 'einaga1979');
    const mh = core.fitMH(anchors);
    const trainEta = core.predictMH(trainSample.Mw_g_mol, mh);
    const testEta = core.predictMH(testSample.Mw_g_mol, mh);
    const rows = trainSample.rows;
    const u = rows.map(r => core.scale(trainEta, r.c_g_L).u);
    const u2 = u.map(v => v ** 2);
    const remainder = rows.map((r, i) => r.eta_rel - 1 - u[i]);
    const kH = dot(u2, remainder) / dot(u2, u2);
    const model = trainResidual(trainSample, rows, trainEta, kH, validateConfig(config));
    return {...model, testSeriesId: testSample.id, strategy: 'wholeGrade', eta_mL_g: testEta,
      trainingEta_mL_g: trainEta, etaPriorSource: 'Einaga 1979 · altı ayrı PS fraksiyonu · Mark–Houwink',
      etaPriorDomain: testSample.Mw_g_mol < mh.minMw || testSample.Mw_g_mol > mh.maxMw ? 'extrapolation' : 'interpolation',
      anchorIds: anchors.map(r => r.row_key), minAnchorMw: mh.minMw, maxAnchorMw: mh.maxMw,
    };
  }
  function predict(model, c_g_L) {
    const {u} = core.scale(model.eta_mL_g, c_g_L);
    const baseline = baselineAt(u, model.kH);
    const cov = model.uTrain.map(v => kernel(u, v, model.settings));
    const latent = forward(model.lower, cov);
    const variance = Math.max(0, kernel(u, u, model.settings) - dot(latent, latent));
    const gp = baseline + dot(cov, model.weights);
    const slope = 1 + 2 * model.kH * u + dot(model.uTrain.map((v, i) => {
      if (u === 0) return 0;
      return cov[i] * (2 / u - (u - v) / model.settings.length_u ** 2);
    }), model.weights);
    const domain = u >= model.minU && u <= model.maxU ? 'interpolation' : 'extrapolation';
    const physicallyAdmissible = Number.isFinite(gp) && gp >= 1 && slope >= 0;
    const eligible = physicallyAdmissible && domain === 'interpolation' && u < 1 && model.etaPriorDomain !== 'extrapolation';
    return {c_g_L, u, baseline, gp, gpConditionalSD: Math.sqrt(variance), slope_dEtaRel_du: slope,
      domain, physicallyAdmissible, eligible, etaPriorDomain: model.etaPriorDomain};
  }
  function metrics(rows, field) {
    let ss = 0, rr = 0, absolute = 0;
    for (const r of rows) {
      const d = r[field] - r.observed;
      ss += d ** 2; rr += (d / r.observed) ** 2; absolute += Math.abs(d);
    }
    return {n: rows.length, mae: absolute / rows.length, rmse: Math.sqrt(ss / rows.length), rmse_pct: 100 * Math.sqrt(rr / rows.length)};
  }
  function score(model, sample, test) {
    const folds = test.map(r => {
      if (!Number.isFinite(r.eta_rel) || r.eta_rel <= 0) throw new RangeError('Değerlendirme için pozitif ölçüm gerekli.');
      return {id: r.id, key: sample.id + '/' + r.id, observed: r.eta_rel, ...predict(model, r.c_g_L)};
    });
    return {seriesId: sample.id, trainSeriesId: model.trainSeriesId, testSeriesId: sample.id, strategy: model.strategy,
      trainIds: model.trainIds, testIds: folds.map(r => r.key), etaPrior_mL_g: model.eta_mL_g, etaPriorSource: model.etaPriorSource,
      etaPriorDomain: model.etaPriorDomain, kH: model.kH, trainRange_u: [model.minU, model.maxU], folds,
      metrics: {baseline: metrics(folds, 'baseline'), gp: metrics(folds, 'gp')},
    };
  }
  function buildBenchmark(finiteSeries, intrinsicRows, config = settings) {
    const within = [];
    for (const sample of finiteSeries) {
      for (const strategy of strategies) {
        const fit = fitWithin(sample, strategy, config);
        within.push(score(fit, sample, splitRows(sample, strategy).test));
      }
    }
    const eligible = finiteSeries.filter(eligibleGrade);
    const transfer = [];
    for (const train of eligible) {
      for (const test of eligible.filter(s => s.id !== train.id && s.Mw_g_mol !== train.Mw_g_mol)) {
        const fit = fitTransfer(train, test, intrinsicRows, config);
        transfer.push(score(fit, test, test.rows));
      }
    }
    return {settings: validateConfig(config), within, transfer, assumptions: {
      target: 'Boyutsuz bağıl viskozite ηᵣ', feature: 'u = c[η]', mean: '1 + u + kH u²',
      kernel: 'σ² u² v² exp(−(u−v)² / (2 ℓ²))',
      conditionalSD: 'Sabit Huggins ortalaması ve sabit çekirdek altında gizli GP düzeltmesinin standart sapması. Ölçüm, mol kütlesi ve model belirsizliği ayrıca değerlendirilir.',
      transfer: 'Aynı çalışma ve çözücü/sıcaklıkta, diğer fraksiyonun yalnız mol kütlesi ve derişimleri kullanılır. Kapiler verilerden türetilen numuneye ait [η] sınama girdisi olarak kullanılmaz.',
      scope: 'Seyreltik, nötr, lineer PS çözeltileri. Derişim ve harici mol kütlesi interpolasyonunda, u < 1 için nicel değerlendirme; mutlak geçerlilik eşiği deneyle belirlenir.',
    }};
  }
  const api = {settings, strategiesTR, fitWithin, fitTransfer, predict, buildBenchmark};
  root.ViscosityGP = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window === 'undefined' ? globalThis : window);
