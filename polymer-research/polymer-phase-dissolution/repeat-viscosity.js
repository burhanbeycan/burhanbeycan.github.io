(() => {
  'use strict';
  const el=id=>document.getElementById(id),catalog=window.THETA_FAMILY_DATA,engine=window.FamilyModels,reader=window.RepeatViscosity;
  if(!catalog?.rows?.length||!engine||!reader||!window.ThetaTransferAI||!window.PolyOcl){el('loadError').hidden=false;return;}
  const ocl=window.PolyOcl(),keys=reader.registry(catalog,ocl),q=new URLSearchParams(location.search);
  const fmt=(n,k=4)=>Number.isFinite(n)?n.toLocaleString('tr-TR',{maximumSignificantDigits:k}):'—';
  const solventNames={'1,4-dioxane':'1,4-dioksan','isoamyl isovalerate':'izoamil izovalerat',cyclohexane:'Siklohekzan'};
  let active=null,model=null,transfer=null;
  for(const [id,v] of Object.entries(catalog.chemistries))el('chemistry').add(new Option(v.label,id));
  function fail(msg){
    active=null;model=null;el('matched').hidden=true;el('issue').textContent=msg;el('issue').hidden=false;
  }
  function draw(){
    const r=catalog.rows.find(v=>v.id===el('sample').value);
    if(!active||!r)return;
    const own=el('massSource').value==='custom',mass=own?1000*el('mass').valueAsNumber:r.mass_g_mol;
    el('massBox').hidden=!own;
    const basis=r.mass_basis==='Mw'?'M𝑤':'Mabs';
    el('massLabel').textContent=basis+' (kg/mol)';
    const t=(r.temperature_basis==='approximate_compilation'?'yaklaşık ':'')+fmt(r.T_C)+' °C';
    const qualifier=r.temperature_basis==='primary'?'birincil yayındaki sıcaklık':'derlemedeki sıcaklık';
    const massText=r.mass_basis==='Mw'?'ağırlık ortalamalı M𝑤':'kaynaktaki mutlak mol kütlesi Mabs; ortalama türü ayrıştırılmamış';
    el('condition').textContent=(solventNames[r.solvent]||r.solvent)+' · '+t+' ('+qualifier+') · kayıtta '+massText+' = '+fmt(r.mass_g_mol/1000,5)+' kg/mol'+(Number.isFinite(r.dispersity)?' · Đ = '+fmt(r.dispersity):'')+'.';
    const sources=el('origin');sources.replaceChildren();
    const source=catalog.sources[r.source];
    const link=document.createElement('a');link.href='https://doi.org/'+source.doi;link.textContent=source.label;sources.append(link);
    const loc=document.createElement('span');loc.textContent=source.location;sources.append(loc);
    if(r.temperature_basis!=='primary'){
      const c=document.createElement('a');c.href=catalog.compilation_url;c.textContent='Deneysel koşulların derlemesi';sources.append(c);
    }
    const err=el('massError');
    if(!(mass>0&&Number.isFinite(mass))){err.textContent='Sonlu, pozitif bir mol kütlesi gir.';err.hidden=false;el('output').hidden=true;return;}
    let out;
    try{out=engine.predict(model,mass,active);}catch(e){err.textContent=e.message;err.hidden=false;el('output').hidden=true;return;}
    err.hidden=true;el('output').hidden=false;
    el('theta').textContent=fmt(out.theta);el('ridge').textContent=fmt(out.ridge);
    el('observed').textContent=own?'—':fmt(r.eta_mL_g);
    el('observedNote').textContent=own?'Kendi numunen için eşlenmiş ölçüm yok.':'Seçilen yayındaki numune · mL/g';
    const target={...r,mass_g_mol:mass};
    // If a custom mass duplicates a catalog sample, exclude that mass too.
    const pool=catalog.rows.filter(ref=>ref.id!==r.id || own);
    const available=own?pool.filter(ref=>Math.abs(ref.mass_g_mol-mass)>1e-8*Math.max(1,mass)):pool;
    if(own)target.id=null;
    const ref=window.ThetaTransferAI.chooseReference(target,available);
    el('refResult').hidden=!ref;
    if(ref){
      const values=window.ThetaTransferAI.predict(transfer,target,ref);
      const src=catalog.sources[ref.source];
      const factor=Math.max(ref.mass_g_mol/mass,mass/ref.mass_g_mol);
      el('refPhysics').textContent=fmt(values.control);el('refAI').textContent=fmt(values.adapted);
      el('refStatus').textContent=ref.sample+' · '+fmt(ref.mass_g_mol/1000,5)+' kg/mol · yayımlanmış [η] = '+fmt(ref.eta_mL_g)+' mL/g. Hedef ve referans aynı nominal kimya, çözücü, sıcaklık ve mol kütlesi tanımını taşır.';
      el('refMethod').textContent='Hedef ile referansın mol kütleleri arasında '+fmt(factor,3)+' kat fark var. İlk değer referansın [η] değerini √(M/Mref) ile ölçekler. İkinci değer, diğer ailelerden öğrenilen ağırlıkla referansın yapı modelinden sapmasını kullanır. Hedefin [η] ölçümü iki hesapta da dışarıda tutulur.';
      el('refSource').textContent=src.label+' · referans ölçüm';el('refSource').href='https://doi.org/'+src.doi;
    } else el('refStatus').textContent='Aynı koşulda kullanılabilecek başka bir referans bulunamadı.';
    const min=model.trainMassMin,max=model.trainMassMax;
    el('support').textContent=engine.labels[model.heldFamily]+' ailesi bütünüyle dışarıda. '+model.nTrain+' eğitim noktasında mol kütlesi '+fmt(min/1000)+'–'+fmt(max/1000)+' kg/mol. Seçilen kütle '+(out.domain==='interpolation'?'aralık içinde':'aralık dışında')+', yapı tanımlayıcısı '+(out.descriptorDomain==='interpolation'?'eğitim aralığında.':'eğitim aralığı dışında.');
  }
  function resolve(sampleId,imported=false){
    let matched;try{matched=reader.match(el('smiles').value,ocl,keys);}catch(e){el('chemistry').value='';fail(e.message);return;}
    active=matched.chemistry;
    const v=catalog.chemistries[active];model=engine.fitForFamily(catalog.rows,v.family);transfer=window.ThetaTransferAI.makeModel(catalog.rows,v.family);
    el('chemistry').value=active;el('name').textContent=v.label;el('issue').hidden=true;el('matched').hidden=false;
    el('structure').innerHTML=reader.depiction(matched);
    el('chemInfo').textContent='Tekrar formülü '+v.repeat_formula+' · '+fmt(v.repeat_molar_mass,6)+' g/mol · tekrar başına '+v.backbone_atoms_per_repeat+' omurga bağı.';
    el('repeatInfo').textContent='Yıldızlar iki bağlantı ucudur. Eşleşen katalog yazımı '+v.repeat_smiles+'.';
    const list=catalog.rows.filter(r=>r.chemistry===active).sort((a,b)=>a.mass_g_mol-b.mass_g_mol);
    if(imported){
      const exactSample=sampleId?list.find(r=>r.id===sampleId):null;
      if(sampleId&&!exactSample){el('chemistry').value='';fail('Aktarılan numune bu tekrar birimiyle eşleşmiyor. Kayıtlı kimyayı seçerek ayrı bir referans hesabı açabilirsin.');return;}
      const reference=exactSample||list[0],inputT=q.get('T');
      const sameT=inputT===null||(inputT.trim()!==''&&Number.isFinite(+inputT)&&Math.abs(+inputT-(reference.T_C+273.15))<1e-5);
      const sameSolvent=!q.has('solvent')||q.get('solvent')===reference.solvent;
      if(!sameT||!sameSolvent){el('chemistry').value='';fail('Aktarılan çözücü veya sıcaklık bu kimyanın kayıtlı θ koşuluyla eşleşmiyor. Kayıtlı bir referans seçebilirsin.');return;}
    }
    const sel=el('sample');sel.replaceChildren();
    for(const r of list)sel.add(new Option(r.sample+' · '+fmt(r.mass_g_mol/1000,5)+' kg/mol · '+catalog.sources[r.source].label,r.id));
    sel.value=list.some(r=>r.id===sampleId)?sampleId:list[Math.floor(list.length/2)].id;
    const row=list.find(r=>r.id===sel.value);el('mass').value=String(row.mass_g_mol/1000);
    el('massSource').value='record';draw();
  }
  el('chemistry').onchange=()=>{const id=el('chemistry').value;el('smiles').value=catalog.chemistries[id]?.repeat_smiles||'';if(id)resolve();else{fail('İki bağlantı noktalı SMILES gir.');el('smiles').focus();}};
  el('smiles').oninput=()=>{el('chemistry').value='';fail('Yeni giriş için Yapıyı eşleştir seçeneğini kullan.');};
  el('smilesForm').onsubmit=e=>{e.preventDefault();resolve();};
  el('sample').onchange=()=>{if(el('massSource').value==='record')el('mass').value=String(catalog.rows.find(r=>r.id===el('sample').value).mass_g_mol/1000);draw();};
  el('massSource').onchange=draw;el('mass').oninput=draw;
  const entry=q.get('repeat'),id=q.get('chemistry'),knownId=Object.hasOwn(catalog.chemistries,id)?id:'PS';
  el('smiles').value=entry!==null?entry:catalog.chemistries[knownId].repeat_smiles;
  el('app').hidden=false;
  if(entry===null&&id!==null&&!Object.hasOwn(catalog.chemistries,id)){el('smiles').value='';fail('Aktarılan kimya θ kohortunda bulunmuyor.');return;}
  if(entry===null&&id===null&&['sample','source','system','Mw','T','solvent'].some(k=>q.has(k))){
    el('smiles').value='';fail('Aktarılan numune için tanınan bir tekrar yapısı gerekli. İki bağlantı noktalı yapıyı gir veya kayıtlı bir kimya seç.');return;
  }
  resolve(q.get('sample'),true);
})();
