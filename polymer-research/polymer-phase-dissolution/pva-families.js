"use strict";
(() => {
  const entries = document.getElementById("entries");
  const status = document.getElementById("status");
  const form = document.getElementById("filters");
  const query = document.getElementById("query");
  const track = document.getElementById("track");
  const evidence = document.getElementById("evidence");
  const catalogue = window.POLYMER_CATALOG;
  if (!Array.isArray(catalogue) || catalogue.length === 0) {
    status.textContent = "Katalog verileri yüklenemedi. CSV dosyası yukarıdaki bağlantıdan açılabilir.";
    for (const control of form.elements) control.disabled = true;
    return;
  }
  const fold = value => String(value || "").toLocaleLowerCase("tr")
    .replaceAll("ı", "i").normalize("NFD").replace(/[\u0300-\u036f]/g, "");
  function el(tag, value, style) {
    const node = document.createElement(tag);
    if (value) node.textContent = value;
    if (style) node.className = style;
    return node;
  }
  function fillFilter(select, values) {
    for (const value of [...new Set(values)].filter(Boolean).sort((a, b) => a.localeCompare(b, "tr"))) {
      select.add(new Option(value, value));
    }
  }
  fillFilter(track, catalogue.map(row => row.track));
  fillFilter(evidence, catalogue.map(row => row.evidence));
  const searchIndex = new Map(catalogue.map(row => [row, fold(Object.values(row).join(" "))]));
  function card(row) {
    const item = el("details", "", "entry");
    const summary = el("summary");
    const title = row.name.includes(row.id) ? row.name : row.id + " · " + row.name;
    const chipStyle = row.evidence === "Araştırma adayı" ? " pending" :
      row.evidence === "İntrinsik viskozite" ? " aux" : "";
    const label = el("span", row.evidence, "measure" + chipStyle);
    const arrow = el("span", "", "toggle");
    arrow.setAttribute("aria-hidden", "true");
    summary.append(el("span", title, "name"), el("span", row.solvent, "solvent"), label, arrow);
    const info = el("div", "", "details");
    const fields = el("dl");
    for (const [term, value] of [
      ["Deneydeki rolü", row.role], ["Kör test grubu", row.group],
      ["Aktarım yolu", row.track], ["Deney koşulları", row.conditions]
    ]) {
      const field = el("div");
      field.append(el("dt", term), el("dd", value || "Belirtilmedi"));
      fields.append(field);
    }
    info.append(fields);
    if (row.note) info.append(el("p", row.note));
    const ref = typeof row.source === "string" ? row.source.trim() : "";
    const doi = /^10\.\d{4,9}\/\S+$/.test(ref);
    if (doi || /^https:\/\//i.test(ref)) {
      const source = el("a", doi ? "Birincil kaynak · " + ref : "Birincil kaynağı aç");
      source.href = doi ? "https://doi.org/" + ref : ref;
      source.target = "_blank";
      source.rel = "noopener noreferrer";
      info.append(source);
    }
    if (row.id.toUpperCase() === "PIB") {
      const data = el("a", "298,15 K NIST ölçümlerini indir (CSV)");
      data.href = "pib-nist-29815.csv";
      data.setAttribute("download", "");
      info.append(el("br"), data);
    }
    item.append(summary, info);
    return item;
  }
  function render() {
    const terms = fold(query.value).trim().split(/\s+/).filter(Boolean);
    const rows = catalogue.filter(row =>
      (!track.value || row.track === track.value) &&
      (!evidence.value || row.evidence === evidence.value) &&
      terms.every(term => searchIndex.get(row).includes(term))
    );
    const fragment = document.createDocumentFragment();
    if (rows.length) {
      for (const row of rows) fragment.append(card(row));
      const count = new Set(rows.map(row => row.group)).size;
      status.textContent = rows.length + " / " + catalogue.length + " kayıt · " + count + " kimyasal grup";
    } else {
      status.textContent = "0 / " + catalogue.length + " kayıt";
      fragment.append(el("p", "Eşleşen kayıt bulunamadı. Aramayı genişletin veya filtreleri temizleyin.", "empty"));
    }
    entries.replaceChildren(fragment);
  }
  query.addEventListener("input", render);
  track.addEventListener("change", render);
  evidence.addEventListener("change", render);
  form.addEventListener("submit", event => event.preventDefault());
  form.addEventListener("reset", () => {
    requestAnimationFrame(() => { render(); query.focus(); });
  });
  render();
})();
