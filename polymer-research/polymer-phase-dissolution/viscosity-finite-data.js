// Capillary measurements transcribed from the named source tables.
(function (root) {
  "use strict";

  const copenhaver = {
    id: "copenhaver-ps-toluene-30",
    label: "PS–toluen · 30 ± 1 °C · Copenhaver vd.",
    polymer: "Polistiren",
    solvent: "Toluen",
    T_K: 303.15,
    temperatureLabel: "30 ± 1 °C",
    temperatureTolerance_K: 1,
    Mn_g_mol: 170000,
    Mw_g_mol: 350000,
    molarMassNote: "Tedarikçi tarafından verilen numune ortalamaları.",
    measurementNote: "Her derişimi ayrı bir ekip ve viskozimetre ölçtü. Tablodaki bağıl viskozite t/t₀ yaklaşımıyla verildi. İkinci ekibin değeri yayında sorgulanmış ve bu seride korunmuştur.",
    source: {
      label: "Copenhaver vd. (2024) · Tablo 1",
      doi: "10.1021/acs.jchemed.4c00439",
      url: "https://doi.org/10.1021/acs.jchemed.4c00439",
      title: "Classic Ubbelöhde Intrinsic Viscosity Laboratory Exercise Made Simple and Fast",
      authors: "Katie Copenhaver, Samuel Pennell, Samyak Jain, Paul S. Russo, Prateek Verma",
      license: "CC BY 4.0",
      licenseUrl: "https://creativecommons.org/licenses/by/4.0/",
    },
    rows: [],
  };

  // Columns: source team, actual c (g/dL), t₀ and its SD, t and its SD (s), published ηrel.
  const teamMeasurements = [
    [1, 0.4295, 110.23, 0.37, 152.97, 1.37, 1.388],
    [2, 0.938, 148.74, 0.43, 246.33, 0.48, 1.656],
    [3, 1.166, 112.44, 2.33, 246.64, 4.03, 2.193],
    [4, 1.722, 109.68, 0.27, 331.73, 0.15, 3.025],
    [5, 1.993, 107.76, 0.61, 391.18, 15.12, 3.630],
  ];
  copenhaver.rows = teamMeasurements.map(([team, c, t0, sd0, t, sd, ratio]) => ({
    id: "team-" + team,
    c_g_L: Number((10 * c).toFixed(3)),
    eta_rel: ratio,
    solutionTime_s: t,
    solventTime_s: t0,
    solutionTimeSD_s: sd,
    solventTimeSD_s: sd0,
    replicates: 3,
    reportedConcern: team === 2,
  }));

  const phywe = {
    id: "phywe-ps-toluene-25",
    label: "PS–toluen · 25,0 °C · PHYWE",
    polymer: "Polistiren",
    solvent: "Toluen",
    T_K: 298.15,
    temperatureLabel: "25,0 °C",
    Mn_g_mol: null,
    Mw_g_mol: null,
    sourceIntrinsic_L_g: 0.0644,
    molarMassNote: "Numunenin Mₙ ve M𝑤 değerleri bildirilmedi. Kılavuz, aynı serinin intrinsik viskozitesinden Mᵥ ≈ 133 500 g/mol hesaplar.",
    measurementNote: "Eğitim amaçlı cihaz kılavuzunda süreler dört ölçümün ortalamasıdır. Saf toluen için t₀ = 76,0 s. Yoğunluklar eşit kabul edilerek bağıl viskozite t/t₀ ile yaklaşık hesaplanmıştır.",
    source: {
      label: "PHYWE P3010601 · Tablo 1",
      doi: null,
      url: "https://spegroup.ru/upload/wikifiles%202/Experiment_guide_P3010601.pdf",
      title: "Determination of the molecular weight of a polymer from intrinsic viscosity measurements",
      authors: "PHYWE",
    },
    rows: [],
  };

  // Columns: c (mg/mL = g/L), mean solution time (s), published ηrel.
  const phyweMeasurements = [
    [4, 97.4, 1.282],
    [8, 122.1, 1.607],
    [12, 150.4, 1.979],
    [16, 181.4, 2.387],
    [20, 219.6, 2.889],
  ];
  phywe.rows = phyweMeasurements.map(([c, t, ratio]) => ({
    id: "c-" + c,
    c_g_L: c,
    eta_rel: ratio,
    solutionTime_s: t,
    solventTime_s: 76.0,
    solutionTimeSD_s: null,
    solventTimeSD_s: null,
    replicates: 4,
    reportedConcern: false,
  }));

  const luhSource = {
    label: "Hwa Luh (1967) · Tablo 8, s. 57",
    title: "The effect of concentration on the viscosity of dilute and moderately concentrated polymer solutions",
    authors: "Hwa Luh",
    url: "https://scholarsmine.mst.edu/masters_theses/6870/",
    doi: null,
  };
  // Original columns: c in g/dL, solution viscosity in cP, and dimensionless viscosity ratio.
  // Both complete high-mass blocks on the continuation of Table 8 are included.
  const luhSamples = [
    { mass: 633000, dispersity: 1.06, intrinsic_L_g: 0.0662, readings: [
      [0.3371, 0.9613, 1.2529], [0.4382, 1.0279, 1.3396],
      [0.5519, 1.1082, 1.4443], [0.6544, 1.1841, 1.5432], [0.9956, 1.4702, 1.9161],
    ]},
    { mass: 1770000, dispersity: 1.08, intrinsic_L_g: 0.1109, readings: [
      [0.1157, 0.8730, 1.1377], [0.1886, 0.9483, 1.2358],
      [0.2877, 1.0579, 1.3787], [0.3870, 1.1800, 1.5380], [0.4704, 1.2907, 1.6821],
    ]},
  ];
  function luhSeries(sample) {
    return {
      id: "luh1967-ps-cyclohexane-" + sample.mass,
      label: "PS–siklohekzan · 34,5 °C · " + sample.mass.toLocaleString("tr-TR") + " g/mol",
      polymer: "Polistiren",
      solvent: "Siklohekzan",
      T_K: 307.65,
      temperatureLabel: "34,5 °C",
      temperatureTolerance_K: 0.02,
      Mn_g_mol: null,
      Mw_g_mol: sample.mass,
      dispersity: sample.dispersity,
      sourceIntrinsic_L_g: sample.intrinsic_L_g,
      molarMassNote: "Pressure Chemical polistiren numunesi. Tablo 1’deki ağırlık ortalamalı mol kütlesi ve dispersite.",
      measurementNote: "Kapiler akış her derişimde beş kez zamanlandı. Dinamik viskozite, akış süresi, kalibrasyon ve yoğunluktan hesaplandı. Tablo 8 hem bu değeri hem bağıl viskoziteyi verir. Tekil süreler ve kayma hızı taraması yayımlanan tabloda yer almıyor.",
      source: luhSource,
      rows: sample.readings.map(([c, eta, ratio], index) => ({
        id: "c-" + (index + 1),
        c_g_L: Number((10 * c).toFixed(3)),
        eta_rel: ratio,
        solutionViscosity_mPa_s: eta,
        solutionTime_s: null,
        solventTime_s: null,
        solutionTimeSD_s: null,
        solventTimeSD_s: null,
        replicates: 5,
        reportedConcern: false,
      })),
    };
  }

  const data = { series: [copenhaver, phywe, ...luhSamples.map(luhSeries)] };
  root.VIS_FINITE = data;
  if (typeof module !== "undefined" && module.exports) module.exports = data;
})(typeof window !== "undefined" ? window : globalThis);
