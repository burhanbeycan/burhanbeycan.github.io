/* Nakamura, Akashi, Norisuye, Teramoto and Sato (1997).
   Only original SAXS results from the unfootnoted rows of Table I, p. 474. */
(function (root) {
  'use strict';
  const source = {
    id: 'nakamura1997', doi: '10.1007/s002890050075',
    url: 'https://doi.org/10.1007/s002890050075',
    caption: 'Nakamura ve ark. (1997) · Tablo I, s. 474',
    method: 'Küçük açılı X-ışını saçılması (SAXS)',
    temperature_C: 34.5, solvent: 'Siklohekzan', temperature_K: 307.65
  };
  const rows = [
    {id:'nakamura1997-F1-B', label:'F1-B', Mw_g_mol:10300, Rg_nm:2.8},
    {id:'nakamura1997-2b-B', label:'2b-B', Mw_g_mol:19800, Rg_nm:4.0}
  ].map(r => ({...r, sourceId:source.id, doi:source.doi, T_K:source.temperature_K, role:'external'}));
  const data = {source, rows};
  root.CHAIN_EXTERNAL = data;
  if (typeof module !== 'undefined' && module.exports) module.exports = data;
})(typeof window === 'undefined' ? globalThis : window);
