(() => {
  const fmt = (n, digits=1) => new Intl.NumberFormat("tr-TR", {
    minimumFractionDigits:digits, maximumFractionDigits:digits
  }).format(n);
  const conc = n => new Intl.NumberFormat("tr-TR",{maximumFractionDigits:3}).format(n);
  const ref = (url,text) => `<a href="${url}" target="_blank" rel="noreferrer">${text}</a>`;

  // The empirical baseline is limited to the same commercial grade and optical series.
  function interpolate(anchors,c) {
    if (!Number.isFinite(c) || c<anchors[0].concentrationGL || c>anchors.at(-1).concentrationGL) return null;
    const j=anchors.findIndex(r=>r.concentrationGL>=c);
    if (j===0 || anchors[j].concentrationGL===c) return anchors[j].tcpC;
    const a=anchors[j-1], b=anchors[j];
    return a.tcpC+(b.tcpC-a.tcpC)*Math.log(c/a.concentrationGL)/Math.log(b.concentrationGL/a.concentrationGL);
  }
  function table(labels,rows) {
    return `<div class="px-scroll"><table><thead><tr>${labels.map(x=>`<th scope="col">${x}</th>`).join("")}</tr></thead><tbody>${rows.map(row=>`<tr>${row.map(x=>`<td>${x}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`;
  }
  function error(anchors,records) {
    const values=records.map(r=>[interpolate(anchors,r.concentrationGL),r.tcpC]).filter(([p])=>p!==null);
    return values.reduce((s,[p,y])=>s+Math.abs(p-y),0)/values.length;
  }
  function chart(records,anchors,external,mode,selected) {
    const x=c=>62+Math.log10(c/0.1)/2*403;
    const y=t=>163-(t-61)*4.35;
    const txt=(x,y,t,anchor="middle")=>`<text x="${x}" y="${y}" text-anchor="${anchor}" fill="#395348" font-family="system-ui,sans-serif" font-size="13">${t}</text>`;
    const guide=(x1,y1,x2,y2)=>`<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="#d1dbd0"/>`;
    const marker=r=>`<circle cx="${x(r.concentrationGL)}" cy="${y(r.tcpC)}" r="4.4" fill="${r.role==="anchor"?"#176456":"white"}" stroke="#176456" stroke-width="2"><title>${conc(r.concentrationGL)} g/L · ${conc(r.tcpC)} °C</title></circle>`;
    const diamond=r=>`<rect x="${x(r.concentrationGL)-3.5}" y="${y(r.tcpC)-3.5}" width="7" height="7" transform="rotate(45 ${x(r.concentrationGL)} ${y(r.tcpC)})" fill="#ad5b28"><title>Zhang · ${conc(r.concentrationGL)} g/L · ${conc(r.tcpC)} °C</title></rect>`;
    const modelValue=interpolate(anchors,selected);
    return `<figure class="px-chart"><svg viewBox="0 0 490 200" aria-label="Sulu Aquazol 50 için derişime göre bulanma sıcaklığı" role="img">
      ${txt(7,15,"T₍cp₎ (°C)","start")}
      ${[65,70,75,80,85,90].map(t=>guide(62,y(t),465,y(t))+txt(48,y(t)+4,String(t),"end")).join("")}
      ${[.1,.5,1,2,5,10].map(c=>guide(x(c),32,x(c),164)+txt(x(c),180,conc(c))).join("")}
      <polyline points="${anchors.map(r=>x(r.concentrationGL)+","+y(r.tcpC)).join(" ")}" fill="none" stroke="#176456" stroke-width="2"/>
      ${records.map(marker).join("")}${mode==="external"?external.map(diamond).join(""):""}
      ${selected!==null && modelValue!==null?`<circle cx="${x(selected)}" cy="${y(modelValue)}" r="8" fill="none" stroke="#27354e" stroke-width="1.4"/>`:""}
      ${txt(465,198,"Derişim (g/L, logaritmik)","end")}</svg>
      <figcaption>Dolu daireler bağlama noktaları, boş daireler aynı seride hesap dışında tutulan ölçümlerdir. Çizgi üç bağlama noktasının logaritmik derişimdeki interpolasyonudur.${mode==="external"?" Turuncu eşkenar dörtgenler Zhang’ın ayrı ısıtma protokolüne aittir.":""}</figcaption></figure>`;
  }
  window.renderPetoxPanel = function(element) {
    const {sources,records,repeat,thermalRecords,externalRecords}=window.PETOX_DATA;
    const [optical,thermal,external]=sources;
    const anchors=records.filter(r=>r.role==="anchor"), tests=records.filter(r=>r.role==="interiorCheck");
    element.classList.add("petox-panel");
    element.innerHTML=`<h2>PEtOx’un suda bulanması</h2><p>Aquazol 50 için deneysel derişim serisini, farklı bir laboratuvarın ölçümlerini ve başka PEtOx numunelerinin kalorimetrik geçişlerini inceleyin.</p>
      <div class="px-tabs" role="group" aria-label="PEtOx veri görünümü">
        <button type="button" data-mode="optical">Bulanma eğrisi</button><button type="button" data-mode="external">Laboratuvarlar arası karşılaştırma</button><button type="button" data-mode="thermal">Kalorimetri</button>
      </div><div class="px-body"></div>`;
    const body=element.querySelector(".px-body");
    function display(mode) {
      for(const b of element.querySelectorAll("[data-mode]")) b.setAttribute("aria-pressed",String(b.dataset.mode===mode));
      if(mode==="optical") {
        body.innerHTML=`<h3>Aquazol 50 · saf su</h3><p>Nominal M<sub>w</sub> 50 kg/mol, Đ 3–4. Hacioglu ve ark. 600 nm’de geçirgenliği izledi. Isıtma hızı 1 °C/dk. T<sub>cp</sub>, geçiş başlangıcı ile sonunun sıcaklık ortalamasıdır.</p>
          <div class="px-graph"></div><div class="px-control"><label>Polimer derişimi (g/L) <input type="number" min="0.1" max="5" step="any" inputmode="decimal" value="1"></label></div>
          <output class="px-value" aria-live="polite"></output><p>Modelin kapsamı bu numunede 0,1–5 g/L arasındaki ısıtma deneyidir. Derişim kütle/hacim cinsindedir.</p>
          <details><summary>Ölçümler ve seri içi kontrol</summary>${table(["Derişim (g/L)","Ölçülen T<sub>cp</sub> (°C)","Hesaplanan (°C)","Kullanım"],records.map(r=>[conc(r.concentrationGL),conc(r.tcpC),fmt(interpolate(anchors,r.concentrationGL)),r.role==="anchor"?"Bağlama":"Kontrol"]))}
          <p>Üç iç kontrol noktasındaki ortalama mutlak fark ${fmt(error(anchors,tests),2)} °C. Bu geriye dönük kontrol aynı numunenin derişim interpolasyonunu sınar. Farklı polimerlere aktarım kapsam dışıdır.</p>
          <p>5 g/L’de kaynak Tablo S2’de 65 °C bildiriyor. Tablo S3’teki ayrı tekrar serisi ${repeat.replicateC.map(x=>fmt(x)).join(" / ")} °C, ortalama ± örnek standart sapması ${fmt(repeat.tcpC)} ± ${fmt(repeat.sampleSdC)} °C (n=3).</p></details>
          <div class="px-links">${ref(optical.url,"Makale · Hacioglu ve ark. (2024)")}${ref(optical.supplement,"Ek tablolar S2 ve S3")}${ref(optical.dataUrl,"Ham veriler · Zenodo")}</div>`;
        const input=body.querySelector("input"),out=body.querySelector("output"),graph=body.querySelector(".px-graph");
        function update() {
          const c=input.valueAsNumber, pred=interpolate(anchors,c);
          graph.innerHTML=chart(records,anchors,externalRecords,mode,Number.isFinite(c)?c:null);
          const obs=records.find(r=>r.concentrationGL===c);
          out.textContent=pred===null?"Bu derişim için hesap aralığı 0,1–5 g/L.":`Hesaplanan ${fmt(pred)} °C${obs?` · Ölçülen ${conc(obs.tcpC)} °C`:" · Ara değer"}`;
        }
        input.addEventListener("input",update); update();
      } else if(mode==="external") {
        body.innerHTML=`<h3>Aynı ticari kalite, iki deney protokolü</h3>
          <p>Hacioglu serisine kurulan eğriyi Zhang’ın hesaplamaya alınmayan Aquazol 50 ölçümleriyle karşılaştırın. Zhang 5 K/dk ısıtmada %50 geçirgenlik sıcaklığını kullanmıştır. Kalite adı aynı olsa da lot kimliği ve uç gruplar eşlenemiyor.</p>
          ${chart(records,anchors,externalRecords,mode,null)}
          ${table(["Derişim (g/L)","Zhang (°C)","Hacioglu modelinden (°C)","Mutlak fark (°C)"],externalRecords.map(r=>{const p=interpolate(anchors,r.concentrationGL);return [conc(r.concentrationGL),fmt(r.tcpC),p===null?"Aralık dışında":fmt(p),p===null?"—":fmt(Math.abs(r.tcpC-p))]}))}
          <p>Aralıkta kalan üç derişimde ortalama mutlak fark ${fmt(error(anchors,externalRecords),2)} °C. Bu değer laboratuvar, numune ve ölçüm tanımındaki değişimleri birlikte içerir. 10 g/L ölçümü hesap aralığını aştığı için dışarıda bırakılmıştır.</p>
          <div class="px-links">${ref(external.supplement,"Zhang (2014) · Tablo 2-2")}${ref(external.relatedArticle,"Zhang ve ark. (2017) · Türbidimetri yöntemleri")}${ref(optical.url,"Hacioglu ve ark. (2024)")}</div>`;
      } else {
        body.innerHTML=`<h3>Serbest PEtOx numunelerinin ısıl geçişi</h3>
          <p>Gauche ve Felisberti, suda 20 g/L serbest PEtOx için yüksek duyarlılıklı diferansiyel taramalı kalorimetri (HSDSC) eğrisinin büküm sıcaklığını raporladı. Kaynak bu sıcaklığı “LCST” olarak etiketler. Ölçülen büyüklük belirli bir derişimdeki kalorimetrik geçiştir.</p>
          ${table(["Numune","M<sub>n</sub> · GPC (kg/mol)","M<sub>w</sub> · GPC (kg/mol)","Đ","HSDSC (°C)"],thermalRecords.map(r=>[r.sample,fmt(r.mnGmol/1000,3),fmt(r.mwGmol/1000,3),fmt(r.dispersity,2),conc(r.transitionC)]))}
          <p>Numuneler 1, 3, 16 ve 24 saatlik sentezlerden geliyor. En uzun iki sentezde GPC’de yüksek mol kütleli omuz ve dallanma bulguları var. Mol kütlesi, zincir mimarisi ve ölçüm yöntemi birlikte değiştiği için bu dört kaydı optik eğriye uydurmak uygun olmaz.</p>
          <div class="px-links">${ref(thermal.url,"Gauche ve Felisberti (2019)")}${ref(thermal.supplement,"Ek Tablo S2 ve Şekil S7")}</div>`;
      }
    }
    for(const b of element.querySelectorAll("[data-mode]")) b.addEventListener("click",()=>display(b.dataset.mode));
    display("optical");
  };
})();
