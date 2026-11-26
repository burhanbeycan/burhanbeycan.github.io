(function (root) {
  "use strict";
  function supported(model, x) {
    if (x.length !== 3 || !x.every(v => Number.isFinite(v) && v > 0)) return false;
    const z = x.map((v, i) => (v - model.reference[i]) / model.step[i]);
    return model.planes.every(([a,b,c,d]) => a*z[0] + b*z[1] + c*z[2] + d <= 1e-8);
  }
  function predict(model, x, key) {
    if (!supported(model, x)) return null;
    if (!(key === "baseline" || key === "quadratic" || key === "ridge")) throw Error("Unknown model.");
    const [u,v,w] = x.map((val,i) => Math.log(val / model.reference[i]));
    const features = [1,u,v,w,u*u,v*v,w*w,u*v,u*w,v*w];
    return Math.exp(model.coef[key].reduce((sum,b,j) => sum + b*features[j], 0));
  }
  // Intersect the vertical (shear-rate) line with the measured design's convex hull.
  function shearInterval(model, composition) {
    const [c,s] = composition.map((v,i) => (v-model.reference[i])/model.step[i]);
    let lo = model.lower[2], hi = model.upper[2];
    for (const [a,b,k,d] of model.planes) {
      const rest = a*c + b*s + d;
      if (Math.abs(k) < 1e-10) { if (rest > 1e-8) return null; continue; }
      const cross = model.reference[2] - model.step[2]*rest/k;
      if (k>0) hi = Math.min(hi,cross); else lo = Math.max(lo,cross);
    }
    if (lo > hi+1e-8) return null;
    return [lo,Math.max(hi,lo)];
  }
  const api={supported,predict,shearInterval};
  root.HPAM=api;
  if (typeof module!=="undefined" && module.exports) module.exports=api;
})(typeof window!=="undefined" ? window : globalThis);
