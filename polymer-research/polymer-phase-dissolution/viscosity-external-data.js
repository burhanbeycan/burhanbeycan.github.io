/* Observations held out by publication and specimen. eta_mL_g is [eta], not eta_rel. */
(function (root) {
  'use strict';
  const turkishStudy = {
    id: 'baysal1992', title: 'Baysal ve Uyanık (1992), Tablolar 1–2',
    doi: '10.1016/0032-3861(92)90695-S', url: 'https://doi.org/10.1016/0032-3861(92)90695-S',
    method: 'Kılcal viskozimetri, c → 0 ekstrapolasyonu', page: 4800
  };
  const standard = {
    id: 'nbs705', title: 'National Bureau of Standards, SRM 705 sertifikası (1978), s. 1',
    url: 'https://tsapps.nist.gov/srmext/certificates/archives/705_COA_1978-11-28.pdf',
    method: 'Sertifikada bildirilen limit viskozite sayısı', page: 1
  };
  const theta = [
    ['PS1', 63000, 1.04, .207], ['PS2', 90000, 1.06, .248], ['PS3', 600000, 1.08, .700],
    ['PS4', 2800000, 1.20, 1.310], ['PS5', 10000000, 1.20, 2.110]
  ].map(([sample, Mw_g_mol, dispersity, eta_dL_g]) => ({
    id:'baysal1992-'+sample, sample, study:turkishStudy.id, system:'ps-cyclohexane-theta',
    solvent:'cyclohexane', T_K:307.65, Mw_g_mol, dispersity, eta_mL_g:100*eta_dL_g,
    training:false, external:true, role:'external', source_url:turkishStudy.url
  }));
  const ref705 = {
    sample:'SRM 705', study:standard.id, Mw_g_mol:179300, Mw_sdm_g_mol:740,
    dispersity:1.07, training:false, external:true, role:'external', source_url:standard.url
  };
  const benzene = {...ref705,id:'nbs705-benzene25',system:'ps-benzene-25',
    solvent:'benzene',T_K:298.15,eta_mL_g:74.3,eta_sdm_mL_g:.18};
  const offTemperature = [
    {...ref705,id:'nbs705-cyclohexane35',system:'ps-cyclohexane-35',solvent:'cyclohexane',
      T_K:308.15,eta_mL_g:35.4,eta_sdm_mL_g:.24,role:'off-temperature'},
    {...ref705,id:'nbs705-benzene35',system:'ps-benzene-35',solvent:'benzene',
      T_K:308.15,eta_mL_g:74.5,eta_sdm_mL_g:.23,role:'off-temperature'}
  ];
  const data = {rows:[...theta,benzene],offTemperature,sources:[turkishStudy,standard]};
  root.VIS_EXTERNAL = data;
  if (typeof module !== 'undefined' && module.exports) module.exports=data;
})(typeof window === 'undefined' ? globalThis : window);
