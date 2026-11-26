/* Solvent-standard-state activity of PS/cyclohexane, Krigbaum & Geymer fraction II.
   Experimental standard-state changes are converted without invoking an EOS. */
(function (root) {
  'use strict';
  const gasConstant = 8.31446261815324; // J mol−1 K−1
  const source = Object.freeze({
    doi: '10.1021/ja01517a022', table: 'III', page: 1864, fraction: 'II',
    Mn_g_mol: 72000, method: 'İzotermal distilasyon', reference: 'Saf sıvı siklohekzan',
    publication: 'W. R. Krigbaum ve D. O. Geymer, JACS 81 (1959), 1859–1868.'
  });
  const records = [
    [34, '0.341', '0.00392'], [34, '0.51', '0.0225'], [34, '0.546', '0.0284'],
    [44, '0.372', '0.0080'], [44, '0.660', '0.087']
  ];
  const rows = Object.freeze(records.map(([C, phiPrinted, minusLogAPrinted]) => Object.freeze({
    id: 'KG-III-F2-' + C + '-' + phiPrinted, T_K: 273.15 + C, phi: Number(phiPrinted),
    minusLogA: Number(minusLogAPrinted), phiPrinted, minusLogAPrinted
  })));
  const reference = root.DENSE_RESULTS;
  const dense = root.DenseFH;
  const sample = reference.sample;
  // Frozen all-eight Table I fit saved as the model used to score the 44 °C pressure series.
  const pressureFit = reference.folds.find(fold => fold.id === 'hold44');
  const coefficients = {baseline: 0, extended: pressureFit.D, guard: pressureFit.D_guard};
  const stateCache = new Map();

  function D(model = 'extended') {
    if (!Object.hasOwn(coefficients, model)) throw new RangeError('Model tanınmıyor.');
    return coefficients[model];
  }
  function convert(T, minusLogA) {
    return {activity: Math.exp(-minusLogA), lnA: -minusLogA, minusLogA,
      deltaMu_J_mol: -gasConstant * T * minusLogA};
  }
  function observed(row) { return convert(row.T_K, row.minusLogA); }
  function predict(T, phi, model = 'baseline') {
    if (!Number.isFinite(T) || T <= 0 || !Number.isFinite(phi) || phi < 0 || phi >= 1) {
      throw new RangeError('Sıcaklık pozitif ve polimer hacim kesri 0 ≤ φ < 1 olmalıdır.');
    }
    const coefficient = D(model);
    const key = T + '/' + model;
    if (!stateCache.has(key)) {
      if (stateCache.size >= 30) stateCache.delete(stateCache.keys().next().value);
      stateCache.set(key, dense.state(T, sample.Mn_g_mol, sample.vbar_ml_g, sample.Vs_ml_mol, reference.base, coefficient));
    }
    const result = dense.evaluate(1000 * phi / sample.vbar_ml_g, stateCache.get(key));
    // Incompressible solvent reference, 1 kPa mL = 0.001 J.
    const minusLogA = result.pressure * sample.Vs_ml_mol / (1000 * gasConstant * T);
    return {...convert(T, minusLogA), Pi_kPa: result.pressure, phases: result.phase.count};
  }
  function compare() {
    return rows.map(row => ({...row, observed: observed(row),
      baseline: predict(row.T_K, row.phi, 'baseline'), extended: predict(row.T_K, row.phi, 'extended'),
      guard: predict(row.T_K, row.phi, 'guard')}));
  }
  function metrics(selected = compare()) {
    function errors(model) {
      if (!selected.length) return {maeMinusLogA: null, maeActivity: null};
      const mae = field => selected.reduce((sum, row) => sum + Math.abs(row[model][field] - row.observed[field]), 0) / selected.length;
      return {maeMinusLogA: mae('minusLogA'), maeActivity: mae('activity')};
    }
    return {n: selected.length, baseline: errors('baseline'), extended: errors('extended'), guard: errors('guard')};
  }

  const api = {source, rows, observed, predict, compare, metrics, D, referenceVolumeMlMol: sample.Vs_ml_mol};
  root.SolventActivity = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window === 'undefined' ? globalThis : window);
