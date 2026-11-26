// Numerical measurements and source metadata; each row is one original sample.
(function (root) {
  'use strict';
  // Table 3: use absolute multiangle light-scattering Mw (the bold column).
  const nagy = [
    ['super-low-88',20100,88,28.7], ['low-88',43600,88,42.6],
    ['medium-low-88',80300,88,65.8], ['medium-88',127000,88,83.3],
    ['high-88',162000,88,101.0], ['super-low-98',23900,98,34.3],
    ['low-98',35800,98,44.3], ['medium-98',101000,98,84.7],
    ['high-98',155000,98,106.9], ['medium-low-92',91900,92,71.1],
    ['high-92',169000,92,106.2], ['medium-96',105000,96,88.2]
  ].map(([grade,Mw_g_mol,DH_pct,eta_mL_g]) => ({grade,Mw_g_mol,DH_pct,eta_mL_g,
    solvent:'aqueous-sodium-nitrate',T_K:308.15,study:'Nagy 2003/2004'}));
  const gebben = [
    ['PVA-96-88',82500,89,79.5,2.06,'25 °C Şekil 1 başlangıç noktasıyla uyumlu. Tablo 1 sıcaklığı ayrı bildirmiyor.'],
    ['PVA-126-98',101000,97.7,105.1,1.60,'25 °C Tablo 2 ile ayrıca doğrulanıyor.']
  ].map(([grade,Mw_g_mol,DH_pct,eta_mL_g,dispersity,conditionNote]) => ({
    grade,Mw_g_mol,DH_pct,eta_mL_g,dispersity,conditionNote,solvent:'water',study:'Gebben 1985'}));
  const korelc = [
    ['PVA88',24000,87,39],['PVA98',23000,98,41]
  ].map(([grade,Mw_g_mol,DH_pct,eta_mL_g]) => ({grade,Mw_g_mol,DH_pct,eta_mL_g,
    solvent:'aqueous-sodium-chloride-azide',T_K:null,study:'Korelc 2024',
    conditionNote:'10 mM NaCl ve %0,002 NaN₃ içeren sulu SEC. SEC sıcaklığı bildirilmedi.'}));
  const urls = {
    parisi:'https://doi.org/10.1122/8.0000435',
    nagy:'https://www.marinalg.org/wp-content/uploads/2012/09/Wu-2004-Handbook-of-size-exclusion-ch10.pdf',
    gebben:'https://doi.org/10.1016/0032-3861(85)90295-2',
    korelc:'https://doi.org/10.1016/j.ejps.2023.106645'
  };
  const data = {nagy,gebben,korelc,urls};
  root.PvohTransferData = data;
  if (typeof module !== 'undefined' && module.exports) module.exports = data;
})(typeof window === 'undefined' ? globalThis : window);
