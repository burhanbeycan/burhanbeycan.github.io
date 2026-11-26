// Incompressible binary Flory–Huggins free energy per kBT per solvent lattice site.
// The integral (free-energy) chi is c0+c1*phi, not the differential/activity chi.
(function (root) {
  "use strict";
  function xlogx(x) { return x === 0 ? 0 : x * Math.log(x); }
  function f(p, N, c0, c1 = 0) {
    return xlogx(p) / N + xlogx(1 - p) + (c0 + c1 * p) * p * (1 - p);
  }
  function df(p, N, c0, c1 = 0) {
    return Math.log(p) / N + 1 / N - Math.log1p(-p) - 1
      + c0 * (1 - 2 * p) + c1 * p * (2 - 3 * p);
  }
  function d2(p, N, c0, c1 = 0) {
    return 1 / (N * p) + 1 / (1 - p) - 2 * c0 + 2 * c1 * (1 - 3 * p);
  }
  function solventPotential(p, N, c0, c1 = 0) {
    // f - phi*f' ; series retains precision in the dilute phase.
    const logarithmic = p < 1e-3 ? -p * p * (0.5 + p * (1 / 3 + p * (0.25 + p / 5))) : Math.log1p(-p) + p;
    return logarithmic - p / N + p * p * (c0 - c1 + 2 * c1 * p);
  }
  // Return the zero of a monotone increasing continuous function in a bracket.
  function bisect(fn, low, high, iterations = 56) {
    for (let i = 0; i < iterations; i++) {
      const m = (low + high) / 2;
      if (m === low || m === high) break;
      if (fn(m) > 0) high = m; else low = m;
    }
    return (low + high) / 2;
  }
  function criticalPhi(N, c1 = 0) {
    if (!c1) return 1 / (1 + Math.sqrt(N));
    return bisect(p => -1 / (N * p * p) + 1 / (1-p) ** 2 - 6*c1, 0, 1);
  }
  function criticalChi(N, c1 = 0) {
    const p = criticalPhi(N, c1);
    return 0.5 * (1 / (N * p) + 1 / (1 - p)) + c1 * (1 - 3 * p);
  }
  // First find a candidate supporting chord on a discrete lower convex envelope.
  // Linear points span the full range, and a graded grid resolves the spinodal pair.
  function hullSlope(N, c0, c1, left, right, count) {
    const sites = [0,1,left,right];
    for (let i=1;i<count;i++) {
      const u=i/count;
      sites.push(u,left*u*u,right+(1-right)*u*u,left+(right-left)*u);
      sites.push(Math.max(0,left-(right-left)*u),Math.min(1,right+(right-left)*u));
    }
    const sorted=[...new Set(sites)].sort((a,b)=>a-b), lower=[];
    for (const p of sorted) {
      const z=[p,f(p,N,c0,c1)];
      while (lower.length>1) {
        const a=lower[lower.length-2], b=lower[lower.length-1];
        if ((b[1]-a[1])*(p-b[0])<(z[1]-b[1])*(b[0]-a[0])) break;
        lower.pop();
      }
      lower.push(z);
    }
    for(let i=1;i<lower.length;i++) {
      const [p0,f0]=lower[i-1], [p1,f1]=lower[i];
      if(p0<=left && p1>=right) return (f1-f0)/(p1-p0);
    }
    return NaN;
  }
  function equilibrium(N, c0, c1 = 0, hullResolution = 32) {
    if (!(Number.isFinite(N) && N > 0 && Number.isFinite(c0) && Number.isFinite(c1))) throw new RangeError("Geçersiz model parametresi");
    const pc = criticalPhi(N, c1), cc = criticalChi(N, c1);
    const base = { N, c0, c1, pc, cc, spin: null, bin: null, diluteUnderflow: false };
    if (c0 <= cc + 1e-13) return base;
    const left = bisect(p => -d2(p,N,c0,c1), 0, pc), right = bisect(p => d2(p,N,c0,c1), pc, 1);
    const slope = p => df(p,N,c0,c1);
    const logL = Math.log(left);
    function phases(m) {
      // Work in log(phi) on the dilute branch: phi can be far below a uniform grid.
      // If needed, use the analytic phi -> 0 limit for log(phi), beyond float range.
      const logEstimate = N * (m + 1 - c0) - 1;
      const logMin = Math.min(-42, logEstimate - 3, logL - 5);
      const lnA = bisect(logP => {
        const p = Math.exp(logP);
        return logP/N + 1/N - Math.log1p(-p) - 1 + c0*(1-2*p) + c1*p*(2-3*p) - m;
      }, logMin, logL, 51);
      const a = Math.exp(lnA), b = bisect(p => slope(p) - m, right, 1, 51);
      return [a, b, lnA];
    }
    const residual = s => {
      const [a,b] = phases(s);
      return solventPotential(a,N,c0,c1) - solventPotential(b,N,c0,c1);
    };
    let lo=slope(right), hi=slope(left);
    const initial=hullSlope(N,c0,c1,left,right,hullResolution);
    if (initial>lo && initial<hi) {
      if (residual(initial)>0) hi=initial; else lo=initial;
    }
    const m=bisect(residual,lo,hi,53);
    const [a,b,lnA] = phases(m);
    return { ...base, spin: [left,right], bin: [a,b], slope: m, logDilute: lnA, diluteUnderflow: a === 0 };
  }
  function classify(p, eq) {
    if (!(p >= 0 && p <= 1)) throw new RangeError("Hacim kesri [0,1] aralığında olmalıdır.");
    const curvature = d2(p,eq.N,eq.c0,eq.c1);
    const out = (label, count, frac = null) => ({ label, count, curvature, richFraction: frac });
    if (!eq.bin) return out(Math.abs(eq.c0 - eq.cc) <= 1e-13 && Math.abs(p-eq.pc)<1e-10 ? "critical" : "one", 1);
    const [a,b] = eq.bin, [l,r] = eq.spin;
    const on = (v,target) => Math.abs(v-target) <= 32*Number.EPSILON*Math.max(1e-15,Math.abs(v),Math.abs(target));
    if (p === 0 || p === 1) return out("one",1);
    if (on(p,a) || on(p,b)) return out("boundary",1);
    if (p < a || p > b) return out("one",1);
    const fraction = (p-a)/(b-a);
    if (on(p,l)||on(p,r)) return out("spinodal",2,fraction);
    return out(p > l && p < r ? "unstable" : "meta",2,fraction);
  }
  const api = {f,df,d2,solventPotential,criticalPhi,criticalChi,equilibrium,classify};
  root.FH = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window === "undefined" ? globalThis : window);
