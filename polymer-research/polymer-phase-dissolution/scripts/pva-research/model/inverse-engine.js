/* Polymer Research Studio • steady-shear PVA/water, Merck 141351 • model version 1.
   Source of training measurements: https://doi.org/10.3390/polym15071743 (Table 2). */
(function (root) {
  "use strict";

  function bracket(values, target) {
    let lower = -Infinity;
    let upper = Infinity;
    for (const value of values) {
      if (value <= target + 1e-8 && value > lower) lower = value;
      if (value >= target - 1e-8 && value < upper) upper = value;
    }
    if (!Number.isFinite(lower) || !Number.isFinite(upper)) return [];
    return lower === upper ? [lower] : [lower, upper];
  }

  function support(windows, wt, temperatureK, rate) {
    const concentrations = bracket(windows.map(row => row[0]), wt);
    const temperatures = bracket(windows.map(row => row[1]), temperatureK);
    if (!concentrations.length || !temperatures.length) {
      return "outside_concentration_or_temperature_range";
    }
    for (const c of concentrations) {
      for (const t of temperatures) {
        const row = windows.find(w => w[0] === c && w[1] === t);
        if (!row) return "incomplete_supporting_grid";
        if (rate < row[2] || rate > row[3]) return "outside_shared_measured_shear_window";
      }
    }
    return "supported_by_neighboring_curve_windows";
  }

  function estimate(state, wt, temperatureK, rate) {
    const c = Math.log(wt / state.reference_concentration_wt_percent);
    const t =
      1000 * (1 / temperatureK - 1 / state.reference_temperature_K) /
      state.inverse_temperature_kernel_scale;
    const row = [1, c, t, c * c, c * t];
    const pars = [0, 0, 0];
    row.forEach((value, i) => {
      for (let j = 0; j < 3; j++) pars[j] += value * state.conditional_parameter_coefficients[i][j];
    });
    const shear = Math.log(rate / state.reference_rate_s);
    const exponent = Math.max(0, Math.min(1, shear < 0 ? pars[1] : pars[2]));
    return Math.exp(pars[0] - exponent * shear);
  }

  function predict(state, wt, temperatureK, rate) {
    if (!state || state.version !== 1 || state.protocol !== "pva-water-merck-141351-table2-v1") {
      throw new TypeError("Unsupported PVA model state");
    }
    if (![wt, temperatureK, rate].every(v => typeof v === "number" && Number.isFinite(v) && v > 0)) {
      return { eta_cP: null, supported: false, reason: "invalid_input" };
    }
    const reason = support(state.shear_windows, wt, temperatureK, rate);
    const supported = reason === "supported_by_neighboring_curve_windows";
    return { eta_cP: supported ? estimate(state, wt, temperatureK, rate) : null, supported, reason };
  }

  const api = Object.freeze({ predict });
  root.PvaInverseModel = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(globalThis);
