"use strict";
(() => {
  const el = id => document.getElementById(id);
  const params = new URLSearchParams(location.search);
  const finiteMass = Number(params.get("finiteMass"));
  const hasSpecies = Number.isFinite(finiteMass) && finiteMass > 0;
  const charge = params.has("charge") ? Number(params.get("charge")) : 0;
  const format = value => value.toLocaleString("tr-TR", { maximumSignificantDigits: 5 });
  const valueOf = id => el(id).value === "" ? NaN : Number(el(id).value);
  if (hasSpecies) {
    el("origin").hidden = false;
    el("origin").textContent = "Aktarılan sonlu yapı · " + (params.get("formula") || "yapı kaydı") + " · " + format(finiteMass) + " g/mol.";
    if (charge === 0) el("basis").add(new Option("Saf ve tek uzunlukta oligomer · aktarılan sonlu molekül", "species"));
  }
  function update() {
    const basis = el("basis").value;
    el("mnrow").hidden = basis !== "mn";
    el("mwrow").hidden = basis !== "mw";
    const mn = basis === "mn" ? valueOf("mn") : basis === "species" ? finiteMass : valueOf("mw") / valueOf("dispersity");
    const kelvin = valueOf("temp") + 273.15, c = valueOf("conc");
    const ok = charge === 0 && Number.isFinite(mn) && mn > 0 && (basis !== "mw" || valueOf("dispersity") >= 1) &&
      Number.isFinite(kelvin) && kelvin > 0 && Number.isFinite(c) && c >= 0;
    el("pi").textContent = "—";
    el("detail").textContent = "";
    el("chainLink").hidden = !ok;
    if (ok) {
      const next=new URLSearchParams({source:"osmotic",system:"custom",Mn:String(mn),T:String(kelvin),c:String(c)});
      if (basis === "mw") next.set("Mw",String(valueOf("mw")));
      if (basis === "species") next.set("Mw",String(finiteMass));
      el("chainLink").href="chain.html?"+next;
    }
    if (!ok) {
      el("note").textContent = charge !== 0 ? "Aktarılan molekül net yük taşıyor. Karşı iyonların katkısı için polielektrolit modeli gerekir." : "Geçerli bir mol kütlesi, sıcaklık ve derişim seçin. Dispersite en az 1 olmalıdır.";
      return;
    }
    const slope = 8.31446261815324 * kelvin / mn;
    el("pi").textContent = format(slope * c) + " kPa";
    el("detail").textContent = "Mₙ = " + format(mn) + " g/mol · Sıfır derişimde Π/c = " + format(slope) + " kPa·L/g.";
    el("note").textContent = basis === "species" ? "Hesap, seçtiğin tek moleküler türe aittir. Sonlu derişimdeki etkileşim düzeltmeleri için ayrıca ölçüm gerekir." : "Bu değer, seçili derişimdeki ideal katkıdır. Gerçek basınçla farkı değerlendirmek için aynı sistemin virial verisi gerekir.";
  }
  document.querySelectorAll("input,select").forEach(element => element.addEventListener("input", update));
  el("basis").addEventListener("change", update);
  update();
})();
