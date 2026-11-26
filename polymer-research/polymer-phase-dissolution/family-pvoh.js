// Parisi et al. (2022), Tables I and III. Intrinsic viscosities: 1 dL/g = 100 mL/g.
(function (root) {
  'use strict';
  const source = {
    doi:'10.1122/8.0000435',
    url:'https://doi.org/10.1122/8.0000435',
    open_url:'https://pure.rug.nl/ws/portalfiles/portal/495697175/8.0000435.pdf',
    label:'Parisi ve ark. (2022), Tablo I ve III, s. 1142–1143',
    study:'parisi2022',
    method:'Cannon–Ubbelohde CUC-75, 25 °C. Her derişimde üç akış süresinin ortalaması; Huggins ve Kraemer doğrularının sıfır derişim kesişimi.',
  };
  // Nominal codes differ from Mw values printed in Table I. Order: Mw and Mn (kg/mol), Đ, approximate DH (%).
  const grades = [
    {grade:'171k-80', mw:171.0, mn:112.3, dispersity:1.52, dh:80},
    {grade:'144k-88', mw:144.5, mn:90.9, dispersity:1.59, dh:88},
    {grade:'49k-75', mw:49.5, mn:28.0, dispersity:1.77, dh:75},
    {grade:'48k-88', mw:48.1, mn:30.0, dispersity:1.61, dh:88},
    {grade:'24k-80', mw:24.5, mn:15.1, dispersity:1.63, dh:80},
  ];
  // First kH from Huggins slope; second kH from Kraemer slope plus 1/2 (Table III).
  const measured = {
    water:[[.78,.88,.59,1.28],[.89,.50,.40,1.12],[.24,2.30,1.80,4.24],[.42,.42,.40,2.40],[.25,.49,.44,3.92]],
    dmso:[[1.70,.36,.36,.59],[1.63,.39,.37,.61],[.58,.45,.40,1.72],[.68,.38,.37,1.47],[.40,.38,.36,2.50]],
  };
  const rows=[];
  for(const [solvent, values] of Object.entries(measured)){
    grades.forEach((g,i)=>{
      const [eta,kH,kHFromKraemer,cstar]=values[i];
      rows.push({
        id:'parisi2022/'+g.grade+'/'+solvent, grade:g.grade, family:'pvoh', solvent,
        study:source.study, source_doi:source.doi, T_K:298.15, Mw_g_mol:g.mw*1000,
        Mn_g_mol:g.mn*1000, dispersity:g.dispersity, DH_pct:g.dh, residual_acetate_mol_pct:100-g.dh,
        eta_mL_g:Math.round(100*eta), eta_dL_g:eta, kH,
        kH_from_Kraemer:kHFromKraemer, cstar_reported_g_dL:cstar, eta_u_mL_g:null, Rg_nm:null,
      });
    });
  }
  const data={source, rows};
  root.FAMILY_PVOH=data;
  if(typeof module!=='undefined'&&module.exports)module.exports=data;
})(typeof window==='undefined'?globalThis:window);
