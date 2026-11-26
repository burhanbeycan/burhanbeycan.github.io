/* Binary incompressible polymer solution, molar-volume reference per solvent molecule. */
(function (root) {
  'use strict';
  const fh = root.FH || (typeof require === 'function' ? require('./phase-engine.js') : null);
  const R = 8.31446261815324; // kPa L mol^-1 K^-1

  function state(T, Mn, polymerMlPerGram, solventMlPerMol, fit) {
    if (![T, Mn, polymerMlPerGram, solventMlPerMol].every(x => Number.isFinite(x) && x > 0)
        || ![fit.A, fit.B, fit.C].every(Number.isFinite)) {
      throw new RangeError('Sıcaklık, molar kütle, hacim ve model katsayılarını kontrol edin.');
    }
    const N = Mn * polymerMlPerGram / solventMlPerMol;
    const chi0 = fit.A + fit.B / T;
    const eq = fh.equilibrium(N, chi0, fit.C);
    const A2 = (0.5 - chi0 + fit.C) * polymerMlPerGram ** 2 / solventMlPerMol;
    return {T, Mn, N, v: polymerMlPerGram, Vs: solventMlPerMol, chi0, C: fit.C, A2, eq};
  }

  function evaluate(c, s) {
    const phi = c * s.v / 1000;
    if (!Number.isFinite(phi) || phi < 0 || phi >= 1) throw new RangeError('Polimer hacim kesri 0 ≤ φ < 1 olmalı.');
    const phase = fh.classify(phi, s.eq);
    // Inside the coexistence interval, the convex envelope has a constant intercept.
    const referencePhi = phase.count === 2 ? s.eq.bin[0] : phi;
    const solventMu = fh.solventPotential(referencePhi, s.N, s.chi0, s.C);
    const pressure = -1000 * R * s.T / s.Vs * solventMu;
    const ideal = R * s.T * c / s.Mn;
    const reduced = c > 0 ? pressure / c : R * s.T / s.Mn;
    return {c, phi, pressure, reduced, ideal, phase};
  }

  const api = {R, state, evaluate};
  root.Osmosis = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window === 'undefined' ? globalThis : window);
