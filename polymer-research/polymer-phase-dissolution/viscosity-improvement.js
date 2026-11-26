(() => {
  'use strict';
  const el = id => document.getElementById(id);
  const set = (id, text) => { el(id).textContent = text; };
  const numberFmt = new Intl.NumberFormat('tr-TR', { maximumFractionDigits: 2 });
  const microFmt = new Intl.NumberFormat('tr-TR', { maximumSignificantDigits: 3 });
  const number = value => Math.abs(value) > 0 && Math.abs(value) < .01 ? microFmt.format(value) : numberFmt.format(value);
  const percentResidual = (row, key) => 100 * (row[key] / row.observed - 1);
  const groups = rows => [...new Set(rows.map(row => row.group))];
  const average = values => values.reduce((sum, value) => sum + value, 0) / values.length;
  const data = window.ViscosityImprovement?.protocols;
  let protocol, rows, activeRow;

  const validProtocol = p => p && typeof p.id === 'string' && typeof p.label === 'string' &&
    typeof p.baseline === 'string' && typeof p.candidate === 'string' && typeof p.unit === 'string' &&
    ['MAPE', 'RMSE'].includes(p.metric) && ['pooled', 'macro'].includes(p.aggregation) &&
    Array.isArray(p.rows) && p.rows.length && p.rows.every(r => r && typeof r.group === 'string' &&
      r.group.length && r.id != null && ['observed', 'baseline', 'candidate'].every(
        key => typeof r[key] === 'number' && Number.isFinite(r[key]) && r[key] > 0));

  if (!Array.isArray(data) || !data.length || !data.every(validProtocol)) {
    el('loadError').hidden = false;
    return;
  }

  function score(subset, key) {
    if (protocol.metric === 'MAPE') return average(subset.map(r => Math.abs(percentResidual(r, key))));
    return Math.sqrt(average(subset.map(r => (r[key] - r.observed) ** 2)));
  }

  function aggregate(subset, key) {
    if (protocol.aggregation === 'pooled') return score(subset, key);
    return average(groups(subset).map(group => score(subset.filter(r => r.group === group), key)));
  }

  function scoreText(scoreValue) {
    return protocol.metric === 'MAPE' ? '%' + number(scoreValue) : number(scoreValue) + ' ' + protocol.unit;
  }

  function options(id, items, wanted) {
    const select = el(id);
    select.replaceChildren(...items.map(([value, label]) => {
      const opt = document.createElement('option');
      opt.value = String(value);
      opt.textContent = label;
      return opt;
    }));
    if (wanted != null && items.some(([v]) => String(v) === String(wanted))) select.value = String(wanted);
  }

  function sourceLink(label, url) {
    if (!label || !url) return null;
    let target;
    try { target = new URL(url, window.location.href); } catch (_) { return null; }
    if (!['http:', 'https:'].includes(target.protocol) &&
      !(target.protocol === 'file:' && window.location.protocol === 'file:' && !/^[\w+.-]+:|^\/\//i.test(url))) return null;
    const link = document.createElement('a');
    link.textContent = label;
    link.href = url;
    return link;
  }

  function rememberSelection() {
    const params = new URLSearchParams(window.location.search);
    params.set('protocol', protocol.id);
    const group = el('group').value;
    if (group) params.set('group', group); else params.delete('group');
    params.set('sample', String(activeRow.id));
    const url = new URL(window.location.href);
    url.search = params.toString();
    try { history.replaceState(null, '', url.href); } catch (_) { /* Standalone file readers can disable history. */ }
  }

  function openProtocol(group, sample) {
    protocol = data.find(p => p.id === el('protocol').value);
    options('group', [['', 'Tüm gruplar'], ...groups(protocol.rows).map(g => [g, g])], group);
    set('protocolDescription', protocol.description || '');
    set('scope', protocol.scope || '');
    el('scope').hidden = !protocol.scope;
    set('baselineName', protocol.baseline);
    set('candidateName', protocol.candidate);
    set('legendBase', protocol.baseline);
    set('legendNew', protocol.candidate);
    set('baseDetail', protocol.baseline);
    set('newDetail', protocol.candidate);
    set('baseGroup', protocol.baseline + ' · ' + protocol.metric);
    set('newGroup', protocol.candidate + ' · ' + protocol.metric);
    set('observedHeader', 'Ölçüm (' + protocol.unit + ')');
    set('baseHeader', protocol.baseline + ' (' + protocol.unit + ')');
    set('newHeader', protocol.candidate + ' (' + protocol.unit + ')');
    set('method', protocol.method || '');
    set('limit', protocol.limit || '');
    el('sources').replaceChildren(...(protocol.links || []).map(l => sourceLink(l.label, l.url)).filter(Boolean));
    openGroup(sample);
  }

  function displayScores() {
    const oldValue = aggregate(rows, 'baseline'), newValue = aggregate(rows, 'candidate');
    const delta = newValue - oldValue;
    set('baselineScore', scoreText(oldValue));
    set('candidateScore', scoreText(newValue));
    const outcome = el('changeLabel').parentElement;
    const tolerance = 1e-12 * Math.max(1, oldValue, newValue);
    const change = Math.abs(delta) <= tolerance ? 'same' : delta > 0 ? 'higher' : 'lower';
    outcome.dataset.change = change;
    set('changeLabel', change === 'lower' ? 'Hata azaldı' : change === 'higher' ? 'Hata arttı' : 'Hata değişmedi');
    set('changeScore', number(change === 'same' ? 0 : Math.abs(delta)) + (protocol.metric === 'MAPE' ? ' yüzde puan' : ' ' + protocol.unit));
    const groupSelected = el('group').value;
    const how = !groupSelected && protocol.aggregation === 'macro' ?
      'Her grubun hatası ayrı hesaplandı ve gruplar eşit ağırlıkla ortalandı.' : 'Her ölçüme eşit ağırlık verildi.';
    const metricText = protocol.metric === 'MAPE' ? 'MAPE, ortalama mutlak yüzde hatadır.' : 'RMSE, ortalama kare hatanın kareköküdür.';
    set('aggregation', rows.length + ' ölçüm · ' + metricText + ' ' + how);
  }

  function textCell(tag, content) {
    const cell = document.createElement(tag);
    cell.textContent = content;
    return cell;
  }

  function linkedRow(id, name, fields, action, active) {
    const tr = document.createElement('tr');
    tr.dataset.record = id;
    tr.classList.toggle('current', active);
    const head = document.createElement('th');
    head.scope = 'row';
    const button = document.createElement('button');
    button.type = 'button';
    button.textContent = name;
    button.addEventListener('click', action);
    head.append(button);
    tr.append(head, ...fields.map(field => textCell('td', field)));
    return tr;
  }

  function openGroup(sample) {
    const group = el('group').value;
    rows = group ? protocol.rows.filter(r => r.group === group) : protocol.rows;
    options('sample', rows.map(r => [r.id, (group ? '' : r.group + ' · ') + (r.label || r.id)]), sample);
    displayScores();
    el('groupRows').replaceChildren(...groups(protocol.rows).map(g => {
      const groupRows = protocol.rows.filter(r => r.group === g);
      return linkedRow(g, g,
        [groupRows.length, scoreText(score(groupRows, 'baseline')), scoreText(score(groupRows, 'candidate'))],
        () => { el('group').value = g; openGroup(); el('group').focus(); }, g === group);
    }));
    el('sampleRows').replaceChildren(...rows.map(r => linkedRow(r.id, r.label || r.id,
      [number(r.observed), number(r.baseline), number(r.candidate)], () => {
        el('sample').value = String(r.id); openSample(); el('sample').focus();
      }, false)));
    openSample();
  }

  function openSample() {
    activeRow = rows.find(r => String(r.id) === el('sample').value) || rows[0];
    const labelValue = v => number(v) + ' ' + protocol.unit;
    const detailValue = key => {
      const delta = percentResidual(activeRow, key);
      return labelValue(activeRow[key]) + ' (' + (delta < 0 ? '−' : delta > 0 ? '+' : '') + '%' + number(Math.abs(delta)) + ')';
    };
    set('condition', [activeRow.group, activeRow.condition].filter(Boolean).join(' · '));
    set('observed', labelValue(activeRow.observed));
    set('baseValue', detailValue('baseline'));
    set('newValue', detailValue('candidate'));
    const reference = sourceLink(activeRow.sourceLabel || 'Ölçümün kaynağı', activeRow.sourceUrl);
    el('sampleSource').replaceChildren(...(reference ? [reference] : []));
    for (const row of el('sampleRows').children) row.classList.toggle('current', row.dataset.record === String(activeRow.id));
    drawResiduals();
    rememberSelection();
  }

  function svg(tag, attrs, text) {
    const elem = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [name, value] of Object.entries(attrs)) elem.setAttribute(name, String(value));
    if (text != null) elem.textContent = text;
    return elem;
  }

  function drawResiduals() {
    const width = 580, height = 320;
    const left = 58, right = width - 13, top = 29, bottom = height - 57;
    const xObserved = rows.map(r => Math.log10(r.observed));
    const xMin = Math.min(...xObserved), xMax = Math.max(...xObserved);
    const xPad = Math.max((xMax - xMin) * .07, .035);
    const lo = xMin - xPad, hi = xMax + xPad;
    const allResiduals = rows.flatMap(r => [percentResidual(r, 'baseline'), percentResidual(r, 'candidate')]);
    const minResidual = Math.min(0, ...allResiduals), maxResidual = Math.max(0, ...allResiduals);
    const crudeStep = Math.max((maxResidual - minResidual) / 4, .1);
    const magnitude = 10 ** Math.floor(Math.log10(crudeStep));
    const step = [1, 2, 2.5, 5, 10].find(c => c >= crudeStep / magnitude) * magnitude;
    const minTick = Math.min(-step, Math.floor(minResidual / step) * step);
    const maxTick = Math.max(step, Math.ceil(maxResidual / step) * step);
    const x = observed => left + (Math.log10(observed) - lo) / (hi - lo) * (right - left);
    const y = residual => bottom - (residual - minTick) / (maxTick - minTick) * (bottom - top);
    const nodes = [
      svg('title', {}, 'Her ölçüm için iki modelin (tahmin / ölçüm − 1) × 100 sapması'),
      svg('text', { x: left, y: 16, 'font-size': 13 }, 'Tahminin sapması (%)'),
      svg('text', { x: (left + right) / 2, y: height - 8, 'text-anchor': 'middle', 'font-size': 13 }, 'Ölçülen viskozite (' + protocol.unit + ')')
    ];
    for (let tick = minTick; tick <= maxTick + step / 100; tick += step) {
      const isZero = Math.abs(tick) < .000001 * step;
      nodes.push(svg('line', { x1: left, x2: right, y1: y(tick), y2: y(tick), stroke: isZero ? '#547164' : '#d8e0d7', 'stroke-dasharray': isZero ? '5 3' : 'none' }));
      nodes.push(svg('text', { x: left - 6, y: y(tick) + 4, 'font-size': 12, 'text-anchor': 'end' }, number(isZero ? 0 : tick)));
    }
    const ticks = xMin === xMax ? [xMin] : [0, .25, .5, .75, 1].map(f => xMin + f * (xMax - xMin));
    for (const t of ticks) {
      const val = 10 ** t, pos = x(val);
      nodes.push(svg('line', { x1: pos, x2: pos, y1: bottom, y2: bottom + 5, stroke: '#718578' }));
      nodes.push(svg('text', { x: pos, y: bottom + 20, 'text-anchor': 'middle', 'font-size': 12 }, microFmt.format(val)));
    }
    const ordered = rows.filter(r => r !== activeRow).concat(activeRow);
    for (const r of ordered) {
      const active = r === activeRow;
      const px = x(r.observed);
      const yBase = y(percentResidual(r, 'baseline')), yNew = y(percentResidual(r, 'candidate'));
      const group = svg('g', { opacity: active ? 1 : .72 });
      if (active) group.append(svg('rect', { x: px - 8, y: Math.min(yBase, yNew) - 8, width: 16, height: Math.abs(yBase - yNew) + 16, rx: 7, fill: '#faefc4', stroke: '#a88125', 'stroke-width': .8 }));
      group.append(svg('title', {}, (r.label || r.id) + ' · ' + r.group + ' · ' + number(r.observed) + ' ' + protocol.unit),
        svg('line', { x1: px, x2: px, y1: yBase, y2: yNew, stroke: '#84968b', 'stroke-width': 1.2 }),
        svg('circle', { cx: px, cy: yBase, r: active ? 5 : 4.1, fill: '#fff', stroke: '#96562c', 'stroke-width': 1.8 }),
        svg('path', { d: 'M' + px + ' ' + (yNew - 4.4) + 'l4.4 4.4-4.4 4.4-4.4-4.4z', fill: '#126754' }));
      nodes.push(group);
    }
    el('residualPlot').replaceChildren(...nodes);
  }

  function downloadCsv() {
    const columns = ['protocol', 'group', 'sample_id', 'sample', 'unit', 'observed', 'baseline_model', 'baseline_prediction', 'candidate_model', 'candidate_prediction', 'baseline_signed_error_pct', 'candidate_signed_error_pct', 'condition', 'source_url'];
    const values = rows.map(r => [protocol.id, r.group, r.id, r.label || r.id, protocol.unit, r.observed, protocol.baseline, r.baseline, protocol.candidate, r.candidate, percentResidual(r, 'baseline'), percentResidual(r, 'candidate'), r.condition || '', r.sourceUrl || '']);
    const csv = [columns, ...values].map(line => line.map(v => {
      if (typeof v === 'number') return String(v);
      const quoted = String(v).replace(/"/g, '""');
      return '"' + (/^[=+@-]/.test(quoted) ? "'" : '') + quoted + '"';
    }).join(',')).join('\r\n');
    const url = URL.createObjectURL(new Blob(['\ufeff', csv], { type: 'text/csv;charset=utf-8' }));
    const link = document.createElement('a');
    link.href = url; link.download = 'viskozite-' + protocol.id + '.csv';
    document.body.append(link); link.click(); link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  const search = new URLSearchParams(window.location.search);
  options('protocol', data.map(p => [p.id, p.label]), search.get('protocol'));
  el('protocol').addEventListener('change', () => openProtocol());
  el('group').addEventListener('change', () => openGroup());
  el('sample').addEventListener('change', openSample);
  el('csv').addEventListener('click', downloadCsv);
  openProtocol(search.get('group'), search.get('sample'));
  el('app').hidden = false;
})();
