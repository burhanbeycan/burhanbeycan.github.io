// Reference-volume Flory–Huggins free energy for the fixed, narrow-PS pressure challenge.
(function (root) {
  'use strict';
  const FH = root.FH || (typeof require === 'function' ? require('./phase-engine.js') : null);
  const Osmosis = root.Osmosis || (typeof require === 'function' ? require('./osmotic-engine.js') : null);
  const maxD = 1 / 12; // Sufficient for f'''' > 0 at every N and every 0 < phi < 1.

  function energy(p, s) {
    return FH.f(p, s.N, s.chi0, s.C) + s.D * p ** 3 * (1 - p);
  }
  function derivative(p, s) {
    return FH.df(p, s.N, s.chi0, s.C) + s.D * p * p * (3 - 4 * p);
  }
  function curvature(p, s) {
    return FH.d2(p, s.N, s.chi0, s.C) + 6 * s.D * p * (1 - 2 * p);
  }
  function solventPotential(p, s) {
    return FH.solventPotential(p, s.N, s.chi0, s.C) + s.D * p ** 3 * (3 * p - 2);
  }
  function bisect(increasing, lo, hi) {
    for (let i = 0; i < 55; i++) {
      const m = (lo + hi) / 2;
      if (m === lo || m === hi) break;
      if (increasing(m) > 0) hi = m;
      else lo = m;
    }
    return (lo + hi) / 2;
  }
  function equilibrium(s) {
    if (s.D === 0) return FH.equilibrium(s.N, s.chi0, s.C);
    // With D <= 1/12 the fourth derivative is positive: the curvature has
    // a unique minimum. Thus either f is convex or exactly one supporting tie line exists.
    const pc = bisect(p => -1 / (s.N * p ** 2) + 1 / (1 - p) ** 2
      - 6 * s.C + 6 * s.D * (1 - 4 * p), 0, 1);
    const base = {pc, spin:null, bin:null};
    if (curvature(pc, s) >= -1e-13) return base;
    const l = bisect(p => -curvature(p, s), 0, pc);
    const r = bisect(p => curvature(p, s), pc, 1);
    const logL = Math.log(l);
    function ends(m) {
      const estimate = s.N * (m + 1 - s.chi0) - 1;
      const lnA = bisect(u => {
        const p = Math.exp(u);
        return (u + 1) / s.N - Math.log1p(-p) - 1 + s.chi0*(1-2*p)
          + s.C*p*(2-3*p) + s.D*p*p*(3-4*p) - m;
      }, Math.min(-42, estimate - 3, logL - 5), logL);
      return [Math.exp(lnA), bisect(p => derivative(p, s) - m, r, 1)];
    }
    // On the two outer convex branches, d(omega_lean-omega_rich)/dm = b-a > 0.
    const m = bisect(slope => {
      const [a,b] = ends(slope);
      return solventPotential(a,s) - solventPotential(b,s);
    }, derivative(r,s), derivative(l,s));
    return {pc, spin:[l,r], bin:ends(m), slope:m};
  }
  function state(T, Mn, vbar, Vs, fit, D) {
    if (![T,Mn,vbar,Vs].every(x => Number.isFinite(x) && x > 0)
      || ![fit?.A,fit?.B,fit?.C,D].every(Number.isFinite) || D > maxD)
      throw new RangeError('Geçersiz termodinamik durum veya katsayı.');
    const chi0 = fit.A + fit.B/T;
    const s = {T,Mn,N:Mn*vbar/Vs,v:vbar,Vs,chi0,C:fit.C,D,
      A2:(0.5-chi0+fit.C)*vbar**2/Vs};
    s.eq = equilibrium(s);
    return s;
  }
  function evaluate(c, s) {
    const p = c*s.v/1000;
    if (!Number.isFinite(p) || p < 0 || p >= 1) throw new RangeError('Geçersiz derişim.');
    const bin = s.eq.bin;
    const two = !!bin && p > bin[0] && p < bin[1];
    const ref = two ? bin[0] : p;
    const pressure = -1000 * Osmosis.R * s.T / s.Vs * solventPotential(ref,s);
    const ideal = Osmosis.R*s.T*c/s.Mn;
    return {c,phi:p,pressure,ideal,reduced:c ? pressure/c : Osmosis.R*s.T/s.Mn,
      phase:{count:two ? 2 : 1},curvature:curvature(p,s)};
  }
  const api = {energy,derivative,curvature,solventPotential,state,evaluate,maxD};
  root.DenseFH = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window === 'undefined' ? globalThis : window);
