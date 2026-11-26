// Paired measurements on identical fractions at two solvent/temperature conditions.
(function (root) {
  'use strict';
  const names = ['cyclohexane', 'benzene'];
  function makePairs(rows) {
    const byId = new Map();
    for (const row of rows.filter(row => row.matched_Rg && row.training)) {
      if (!names.includes(row.solvent) || ![row.Rg_nm, row.eta_mL_g, row.Mw_g_mol].every(n => Number.isFinite(n) && n > 0))
        throw new Error('Eksik veya geçersiz eşleşme ölçümü.');
      if (!byId.has(row.id)) byId.set(row.id, {id:row.id, parent:row.parent});
      const pair = byId.get(row.id);
      if (pair[row.solvent] || !row.parent || pair.parent !== row.parent) throw new Error('Fraksiyon kimliği tekil olmalı.');
      pair[row.solvent] = row;
    }
    const pairs = Array.from(byId.values());
    for (const p of pairs) {
      if (!p.benzene || !p.cyclohexane || p.benzene.Mw_g_mol !== p.cyclohexane.Mw_g_mol)
        throw new Error('Aynı fraksiyonun iki koşulu ve ortak nominal mol kütlesi gerekir.');
    }
    return pairs.sort((a,b) => a.benzene.Mw_g_mol - b.benzene.Mw_g_mol);
  }
  function ratio(pair, from, to) {
    return (pair[to].eta_mL_g / pair[from].eta_mL_g) / (pair[to].Rg_nm / pair[from].Rg_nm) ** 3;
  }
  function benchmark(pairs, from, split) {
    if (!names.includes(from) || !['fraction','parent'].includes(split)) throw new Error('Koşul veya bölme geçersiz.');
    const to = names.find(n => n !== from);
    return pairs.map(pair => {
      const train = pairs.filter(p => split === 'fraction' ? p.id !== pair.id : p.parent !== pair.parent);
      if (!train.length) throw new Error('Düzeltme için bağımsız eğitim fraksiyonu gerekir.');
      const factor = Math.exp(train.reduce((v,p) => v + Math.log(ratio(p, from, to)),0) / train.length);
      const source = pair[from], target = pair[to];
      const radiusRatio = target.Rg_nm / source.Rg_nm;
      const unscaled = source.eta_mL_g * radiusRatio ** 3;
      const adjusted = factor * unscaled;
      const observed = target.eta_mL_g;
      return {id:pair.id, parent:pair.parent, from, to, radiusRatio, factor, observed, unscaled, adjusted,
        sourceEta:source.eta_mL_g, sourceRg:source.Rg_nm, targetRg:target.Rg_nm,
        baselineError:100*(unscaled/observed-1), adjustedError:100*(adjusted/observed-1), train:train.map(p=>p.id)};
    });
  }
  const api = {makePairs, benchmark};
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  if (root) root.ViscosityTransfer = api;
})(typeof window === 'undefined' ? null : window);
