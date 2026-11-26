(() => {
  'use strict';
  const el=id=>document.getElementById(id),engine=window.PvohTransfer,source=window.PvohTransferData,rows=window.FAMILY_PVOH?.rows;
  if(!engine||!source||!rows){el('loadError').hidden=false;return;}
  const {folds,external}=engine.evaluate(rows,source),full=engine.fit(rows,source.nagy);
  const names={mass:'Mol kütlesi',sharedDH:'Mol kütlesi + ortak hidroliz etkisi',partialDH:'Çözücüye özgü hidroliz eğimi · birincil aday',
    partialBoth:'Çözücüye özgü kütle ve hidroliz eğimleri',sourceTransfer:'Sulu SEC verisinden eğim aktarımı'};
  const num=(n,d=2)=>Number(n).toLocaleString('tr-TR',{maximumFractionDigits:d});
  const tableRow=(table,values)=>{const tr=table.insertRow();for(const x of values)tr.insertCell().textContent=x;return tr;};
  for(const k of engine.methods)el('method').add(new Option(names[k],k)); el('method').value='partialBoth';
  const group=()=>el('experiment').value,method=()=>el('method').value;
  function records() {
    if(group()!=='paired')return external[group()].records;
    return folds.filter(r=>el('solvent').value==='all'||el('solvent').value===r.solvent);
  }
  function query() {
    const mass=el('queryMass').valueAsNumber,dh=el('queryDH').valueAsNumber,s=el('querySolvent').value;
    if(!Number.isFinite(mass)||!Number.isFinite(dh)||mass<=0||dh<=0||dh>100){
      el('prediction').textContent='Pozitif M𝑤 ve 0–100 aralığında pozitif hidroliz derecesi gerekiyor.';return;
    }
    const p=engine.predict(full,{solvent:s,Mw_g_mol:mass*1000,DH_pct:dh});
    let domain=[];
    for(const [key,value,label] of [['Mw_g_mol',mass*1000,'M𝑤'],['DH_pct',dh,'hidroliz derecesi']]){
      const [lo,hi]=full.trainRange[key];if(value<lo||value>hi)domain.push(label);
    }
    el('prediction').textContent=names[method()]+' için [η] = '+num(p[method()])+' mL/g. '+
      (domain.length?'Eğitim aralığını aşan girdiler: '+domain.join(', ')+'.':'Her iki girdi de eğitim değerlerinin alt ve üst sınırları içinde.');
  }
  function render() {
    el('solvent').disabled=group()!=='paired';if(group()!=='paired')el('solvent').value='water';
    const data=records(), scores=engine.summarize(data), methodKey=method();
    const notes={
      paired:'Parisi ve arkadaşlarının beş PVOH numunesinin her biri, sudaki ve DMSO’daki ölçümleriyle birlikte ayrıldı. Her tahmin dört başka numuneyle eğitilen modelden geldi.',
      gebben:'Gebben ve arkadaşlarının çapraz bağlama öncesindeki iki başlangıç numunesi. M𝑤 ışık saçılması ile belirlendi. İki bağımsız numunenin hidroliz derecesi (%89 ve %97,7) Parisi eğitim aralığının (%75–88) üzerinde.',
      korelc:'Korelc ve arkadaşlarının iki fraksiyonsuz ticari numunesi. 10 mM NaCl ve %0,002 NaN₃ içeren sulu SEC elüentindeki [η], 25 °C su modeliyle karşılaştırıldı. SEC sıcaklığı bildirilmedi. Ortam, yöntem ve laboratuvar birlikte değişiyor.'
    };
    el('protocol').textContent=notes[group()];
    const summary=el('scores'), samples=el('samples');summary.replaceChildren();samples.replaceChildren();
    for(const k of engine.methods){const s=scores[k],tr=tableRow(summary,[names[k],num(s.rmse),num(s.mae),num(s.mape)]);if(k===methodKey)tr.className='highlight';}
    for(const r of data){
      const solvent=r.solvent==='dmso'?'DMSO':r.solvent==='water'?'Su':'Sulu SEC';
      tableRow(samples,[r.grade,solvent,num(r.Mw_g_mol/1000),num(r.DH_pct),num(r.measured),num(r[methodKey]),num(r[methodKey]-r.measured)]);
    }
    const best=scores.mass.rmse<scores.sharedDH.rmse?'mass':'sharedDH';
    el('comparison').textContent='Bu seçimin daha düşük hatalı eski modeli: '+names[best]+' ('+num(scores[best].rmse)+' mL/g). '+names[methodKey]+' için RMSE '+num(scores[methodKey].rmse)+' mL/g.';
    el('scopeNote').textContent=group()==='paired'?'Kütle aralığı 24,5–171 kg/mol. DH aralığı %75–88. İki eğim eklemek DMSO’da hatayı düşürürken suda ortak hidroliz modeli daha düşük hata verir.':group()==='gebben'?'PVA-126-98 için 105,1 mL/g değeri 25 °C’de ayrıca raporlandı. PVA-96-88’in Tablo 1 değeri Şekil 1’deki 25 °C başlangıç noktasıyla uyumlu.':'Korelc çalışmasındaki PVA35, PVA50 ve PVA74 için ölçülmüş [η] bulunmuyor. 10 % çözeltinin Brookfield viskozitesi ayrı bir fiziksel büyüklüktür.';
    query();
  }
  el('download').onclick=()=>{
    const list=[['grade','solvent','Mw_g_mol','DH_pct','observed_mL_g',...engine.methods.map(k=>k+'_mL_g'),'study','condition','training_grades']];
    for(const r of records())list.push([r.grade,r.solvent,r.Mw_g_mol,r.DH_pct,r.measured,...engine.methods.map(k=>r[k]),r.study||'Parisi 2022',
      r.conditionNote||'25 °C',r.trainingGrades?.join('|')||full.trainingGrades.join('|')]);
    const escape=v=>'"'+String(v).replaceAll('"','""')+'"';
    const a=document.createElement('a'),url=URL.createObjectURL(new Blob([list.map(row=>row.map(escape).join(',')).join('\n')+'\n'],{type:'text/csv;charset=utf-8'}));
    a.href=url;a.download='PVOH_'+group()+'_model_sinamasi.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  el('experiment').onchange=()=>{el('solvent').value='all';render();};el('solvent').onchange=render;el('method').onchange=render;
  el('queryMass').oninput=query;el('queryDH').oninput=query;el('querySolvent').onchange=query;
  el('app').hidden=false;render();
})();
