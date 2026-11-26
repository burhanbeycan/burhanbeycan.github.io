(() => {
  'use strict';
  const hydro = window.HYDRO_REFERENCE;
  const el = id => document.getElementById(id);
  const set = (id, text) => { el(id).textContent = text; };
  const fmt = new Intl.NumberFormat('tr-TR', { maximumFractionDigits: 2 });
  const precise = new Intl.NumberFormat('tr-TR', { maximumSignificantDigits: 4 });
  const massKg = value => precise.format(value / 1000);
  const interval = values => values.map(massKg).join('–') + ' kg/mol';
  const svgNS = 'http://www.w3.org/2000/svg';
  if (!hydro || !hydro.targets.length || !Object.keys(hydro.models || {}).length) {
    el('loadIssue').hidden = false;
    return;
  }
  const author = {
    '10.1002/polb.1989.090271209': 'Roovers ve Martin',
    '10.1021/ma00182a043': 'Tsunashima ve ark.',
    '10.1021/ma00176a048': 'Davidson ve ark.',
    '10.1002/pol.1972.160100103': 'Matsumoto ve ark.',
    '10.1021/ma00029a020': 'Nakamura ve ark.',
    '10.1021/ma00175a030': 'Colby ve ark.'
  };
  const opt = (key, label) => {
    const o = document.createElement('option');
    o.value = key; o.textContent = label; return o;
  };
  let chemistry, targets, customBounds;
  const selectedKnown = () => targets.find(t => t.id === el('target').value);
  function changeFamily() {
    chemistry = el('polymer').value;
    targets = hydro.targets.filter(t => t.family === chemistry).sort((a, b) => a.M_g_mol - b.M_g_mol);
    const samples = hydro.catalogue.filter(t => t.family === chemistry);
    const model = hydro.models[chemistry];
    customBounds = [
      Math.max(model.M_min, Math.min(...samples.map(r => r.M_g_mol))),
      Math.min(model.M_max, Math.max(...samples.map(r => r.M_g_mol)))
    ];
    el('target').replaceChildren(...targets.map(t => opt(t.id, t.id + ' · ' + massKg(t.M_g_mol) + ' kg/mol')));
    el('target').value = targets[Math.floor(targets.length / 2)].id;
    el('mass').min = customBounds[0] / 1000;
    el('mass').max = customBounds[1] / 1000;
    el('mass').value = (Math.sqrt(customBounds[0] * customBounds[1]) / 1000).toPrecision(4);
    changeMode();
  }

  function changeMode() {
    const known = el('mode').value === 'record';
    el('knownField').hidden = !known;
    el('customField').hidden = known;
    const basis = chemistry === 'pib' ? 'Ağırlık ortalamalı mol kütlesi (M𝑤)' : 'Derlemedeki mol kütlesi (M)';
    set('massHint', known ? basis + ' gösterilir. Seçilen numunenin tüm yayını referans seçiminden çıkarılır.' :
      basis + ' kullanılır. Bu kimyada referanslarla üs eğitiminin ortak kütle aralığı ' + interval(customBounds) + '.');
    update();
  }

  function failed(message) {
    el('prediction').hidden = true;
    el('entryIssue').textContent = message;
    el('entryIssue').hidden = false;
  }

  function linkedSource(id, prefix, doi) {
    const p = el(id);
    p.replaceChildren();
    if (!doi) return;
    const a = document.createElement('a');
    a.href = 'https://doi.org/' + doi;
    a.textContent = prefix + ' · ' + (author[doi] || doi);
    p.append(a);
  }

  function update() {
    const known = el('mode').value === 'record';
    const target = known ? selectedKnown() : null;
    const mass = target ? target.M_g_mol : el('mass').valueAsNumber * 1000;
    el('entryIssue').hidden = true;
    if (!Number.isFinite(mass) || mass <= 0) return failed('Pozitif bir mol kütlesi gir (kg/mol).');
    if (!known && (mass < customBounds[0] || mass > customBounds[1])) {
      return failed('Girilen mol kütlesi desteklenen ' + interval(customBounds) + ' aralığının dışında.');
    }
    let result;
    try { result = hydro.predict(chemistry, mass, target?.id); }
    catch (e) { return failed(e.message); }
    if (!result.anchor || !(result.calibrated > 0)) {
      return failed('Aynı kütledeki yayın dışarıda bırakıldığında uygun bir bağımsız referans kalmıyor. Başka bir mol kütlesi veya sınanmış bir numune seç.');
    }
    const anchor = result.anchor;
    const signed = target ? 100 * (result.calibrated / target.eta_mL_g - 1) : null;
    set('etaHat', fmt.format(result.calibrated));
    set('comparisonHeader', target ? 'Ölçülen [η]' : 'Hedef mol kütlesi');
    set('comparisonValue', target ? fmt.format(target.eta_mL_g) : massKg(mass));
    set('comparisonNote', target ? 'mL/g · tahminin sapması ' + (signed > 0 ? '+' : signed < 0 ? '−' : '') + '%' + fmt.format(Math.abs(signed)) : 'kg/mol · kendi numunen');
    set('exponent', precise.format(result.slope));
    set('training', Object.entries(hydro.names).filter(([f]) => f !== chemistry).map(([, name]) => name).join(' ve ') + ' ile öğrenildi.');
    set('anchorId', anchor.id); set('anchorM', massKg(anchor.M_g_mol) + ' kg/mol');
    set('anchorEta', fmt.format(anchor.eta_mL_g) + ' mL/g'); set('ratio', precise.format(mass / anchor.M_g_mol) + ' ×');
    linkedSource('anchorSource', 'Referans yayını', anchor.source_doi);
    linkedSource('targetSource', 'Saklı ölçümün yayını', target?.source_doi);
    const asRange = result.referenceRange[0] === result.referenceRange[1] ?
      'Uygun referans tek kütlede (' + massKg(result.referenceRange[0]) + ' kg/mol).' :
      'Uygun referanslar ' + interval(result.referenceRange) + ' aralığında.';
    const flags = [];
    if (result.outsideReferenceRange) flags.push('Hedef kütle uygun referansların aralığı dışında.');
    if (result.outsideTrainingRange) flags.push('Hedef kütle üssün öğrenildiği aralığın dışında.');
    if (!known && result.excludedDois.length) flags.push('Katalogdaki eş kütleli numunelerin yayınları referans seçiminden çıkarıldı.');
    set('domainInfo', asRange + ' Üssün eğitim aralığı ' + interval(result.trainingRange) + '. ' + flags.join(' '));
    el('observedKey').hidden = !target;
    drawCurve(mass, result, target);
    el('prediction').hidden = false;
  }

  function mark(tag, attrs, text) {
    const e = document.createElementNS(svgNS, tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
    if (text != null) e.textContent = text;
    return e;
  }

  function drawCurve(mass, result, target) {
    const a = result.anchor;
    const lo = Math.min(Math.log10(mass), Math.log10(a.M_g_mol));
    const hi = Math.max(Math.log10(mass), Math.log10(a.M_g_mol));
    const xPad = Math.max(.045, (hi - lo) * .08);
    const xmin = lo - xPad, xmax = hi + xPad;
    const curve = m => a.eta_mL_g * (m / a.M_g_mol) ** result.slope;
    const etaValues = [curve(10 ** xmin), curve(10 ** xmax), result.calibrated, a.eta_mL_g];
    if (target) etaValues.push(target.eta_mL_g);
    const min = Math.log10(Math.min(...etaValues)), max = Math.log10(Math.max(...etaValues));
    const yPad = Math.max(.025, (max - min) * .08);
    const ymin = min - yPad, ymax = max + yPad;
    const x = m => 63 + (Math.log10(m) - xmin) / (xmax - xmin) * 463;
    const y = eta => 255 - (Math.log10(eta) - ymin) / (ymax - ymin) * 226;
    const els = [
      mark('title', {}, 'Referans ölçümüne bağlanmış Mark–Houwink kuvvet yasası'),
      mark('text', { x: 63, y: 18, 'font-size': 13 }, 'İntrinsik viskozite (mL/g)'),
      mark('text', { x: 289, y: 308, 'text-anchor': 'middle', 'font-size': 13 }, 'Mol kütlesi (kg/mol)')
    ];
    for (let i = 0; i <= 4; i++) {
      const m = 10 ** (xmin + i * (xmax - xmin) / 4), eta = 10 ** (ymin + i * (ymax - ymin) / 4);
      const px = x(m), py = y(eta);
      els.push(mark('line', { x1: 63, x2: 526, y1: py, y2: py, stroke: '#d8e0d7' }),
        mark('text', { x: 56, y: py + 4, 'font-size': 12, 'text-anchor': 'end' }, fmt.format(eta)),
        mark('line', { x1: px, x2: px, y1: 255, y2: 260, stroke: '#809183' }),
        mark('text', { x: px, y: 278, 'text-anchor': 'middle', 'font-size': 12 }, massKg(m)));
    }
    els.push(mark('line', { x1: 63, x2: 526, y1: y(curve(10 ** xmin)), y2: y(curve(10 ** xmax)), stroke: '#247255', 'stroke-width': 2 }),
      mark('rect', { x: x(a.M_g_mol) - 5, y: y(a.eta_mL_g) - 5, width: 10, height: 10, fill: '#4b637a' }));
    const px = x(mass), py = y(result.calibrated);
    if (target) {
      const pt = y(target.eta_mL_g);
      els.push(mark('line', { x1: px, x2: px, y1: pt, y2: py, stroke: '#846e55', 'stroke-dasharray': '3 3' }),
        mark('circle', { cx: px, cy: pt, r: 7, fill: '#fff', stroke: '#96562c', 'stroke-width': 2 }));
    }
    els.push(mark('path', { d: 'M' + px + ' ' + (py - 5.5) + 'l5.5 5.5-5.5 5.5-5.5-5.5z', fill: '#126754' }));
    el('massDiagram').replaceChildren(...els);
  }

  el('polymer').replaceChildren(...Object.entries(hydro.names).map(([k, name]) => opt(k, name)));
  const family = new URLSearchParams(window.location.search).get('family');
  if (Object.hasOwn(hydro.names, family)) el('polymer').value = family;
  el('polymer').addEventListener('change', changeFamily);
  el('mode').addEventListener('change', changeMode);
  el('target').addEventListener('change', update);
  el('mass').addEventListener('input', update);
  changeFamily();
  el('refApp').hidden = false;
})();
