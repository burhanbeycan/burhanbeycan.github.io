"use strict";
(() => {
  const $ = (id) => document.getElementById(id);
  const known = window.polymerReference;
  const graphSchema = "polymer-studio/browser-graph/1";
  const identitySchema = "polymer-identity/1";
  const nativeSchema = "polymer-structure-energy/1";
  let saved, geo, yaw = 0.45, tilt = 0.32, activeVersion = 0, sketch;
  const normal = value => typeof value === "string" && value.length > 0;
  function requireThat(test, message) { if (!test) throw new Error(message); }
  function status(message = "", error = false) {
    $("status").textContent = message;
    $("status").classList.toggle("error", error);
  }
  function invalidate() {
    activeVersion += 1;
    saved = undefined; geo = undefined; $("output").hidden = true;
    status();
  }
  function tabs(key) {
    document.querySelectorAll("[data-pane]").forEach(p => p.hidden = p.dataset.pane !== key);
    document.querySelectorAll("[data-tab]").forEach(b => b.setAttribute("aria-pressed", String(b.dataset.tab === key)));
    if (key === "sketch" && !sketch) sketch = window.PolySketch.mount($("sketch"));
  }
  document.querySelectorAll("[data-tab]").forEach(b => b.onclick = () => { invalidate(); tabs(b.dataset.tab); });
  $("smiles").oninput = invalidate;
  window.addEventListener("polymer-graph-change", () => { if (!$("sketch").closest("[data-pane]").hidden) invalidate(); });
  const stable = value => {
    if (Array.isArray(value)) return value.map(stable);
    if (value && typeof value === "object") return Object.fromEntries(Object.keys(value).sort().map(k => [k, stable(value[k])]));
    return value;
  };
  async function digest(input) {
    requireThat(globalThis.crypto?.subtle, "Kimlik doğrulaması için HTTPS veya localhost üzerinden açın.");
    const bytes = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(input));
    return Array.from(new Uint8Array(bytes), x => x.toString(16).padStart(2, "0")).join("");
  }
  async function nativeId(label, value) {
    const json = JSON.stringify(stable(value)).replace(/[^\x00-\x7f]/g, chr => "\\u" + chr.charCodeAt(0).toString(16).padStart(4, "0"));
    return label + "-" + await digest(identitySchema + "\n" + json);
  }
  function usableGeometry(g) {
    if (g == null) return null;
    requireThat(g.unit === "angstrom" && Array.isArray(g.atoms) && g.atoms.length > 0 && g.atoms.length <= 1100 && Array.isArray(g.bonds),
      "Koordinat birimi veya atom/bağ kaydı geçersiz.");
    requireThat(g.atoms.every(a => /^[A-Z][a-z]?$/.test(a.element) && [a.x,a.y,a.z].every(Number.isFinite)) &&
      g.bonds.every(b => [b.a,b.b].every(i => Number.isInteger(i) && i >= 0 && i < g.atoms.length) && b.a !== b.b && b.order > 0),
      "Konformer atomları, bağ indeksleri veya koordinatları geçersiz.");
    return g;
  }
  async function acceptNative(data) {
    const structure = data.schema === identitySchema ? data : data.structure;
    requireThat(structure && structure.schema === identitySchema && normal(structure.canonical_smiles) &&
      structure.recipe && ["repeat","sequence","star","finite"].includes(structure.mode) && structure.mode === structure.recipe.mode &&
      structure.tacticity === structure.recipe.tacticity && normal(structure.formula) && Number.isFinite(structure.molar_mass) &&
      Number.isInteger(structure.heavy_atoms) && structure.heavy_atoms > 0 && typeof structure.stereo_resolved === "boolean",
      "Yapı ve Konformer Enerjisi kimlik alanları eksik veya geçersiz.");
    const ids = await Promise.all([nativeId("structure", structure.canonical_smiles), nativeId("recipe", structure.recipe)]);
    requireThat(ids[0] === structure.structure_id && ids[1] === structure.recipe_id,
      "Yapı veya reçete kimliği dosyadaki SMILES ve reçete ile eşleşmiyor.");
    if (data.schema === nativeSchema) {
      requireThat(/^[a-f0-9]{64}$/.test(data.artifact_id) && data.energy?.canonical_smiles === structure.canonical_smiles,
        "Konformer kaydı ile kimyasal yapı eşleşmiyor.");
    }
    const g = usableGeometry(data.energy?.lowest_geometry);
    if (g) requireThat(g.atoms.filter(a => a.element !== "H").length === structure.heavy_atoms,
      "Konformerdeki ağır atom sayısı yapı kaydıyla uyuşmuyor.");
    return {raw:data, structure, g, exact:structure.structure_id === known.structure.structure_id, native:true};
  }
  function parseFinite(smiles) {
    requireThat(normal(smiles) && smiles.length <= 1800 && !/\s/.test(smiles), "Boşluksuz tek bir sonlu SMILES girin (en çok 1800 karakter).");
    requireThat(!smiles.includes("*") && !smiles.includes("{") && !/:[0-9]+\]/.test(smiles), "Bağlantı noktaları veya stokastik birimler için aşağıdaki tekrar birimi uygulamasını kullanın. Sonlu zincirde uçları açıkça tanımlayın.");
    const O = window.PolyOcl();
    let mol;
    try { mol = O.Molecule.fromSmiles(smiles); }
    catch { throw new Error("SMILES okunamadı. Parantezleri, halka kapanışlarını, atomları ve bağları kontrol edin."); }
    const count = mol.getAllAtoms();
    requireThat(count > 0 && mol.getFragments().length === 1, "Tek bir bağlı, sonlu atom grafiği girin.");
    for (let i = 0; i < count; i++) {
      requireThat(mol.getAtomicNo(i) > 0 && mol.getOccupiedValence(i) <= mol.getMaxValence(i), "Atom " + (i + 1) + " izin verilen değerliği aşıyor.");
    }
    try { mol.validate(); }
    catch (err) {
      // OpenChemLib's overall neutrality check is separate from valence validation.
      if (!String(err).includes("unbalanced atom charge")) throw new Error("Kimyasal doğrulama başarısız. Bağ derecelerini ve stereokimyayı kontrol edin.");
    }
    return {O, mol, count};
  }
  async function buildGraph(smiles) {
    const {O, mol, count} = parseFinite(smiles);
    const nativeGraph = mol.getIDCode(), formula = mol.getMolecularFormula();
    const structure = {canonical_smiles:mol.toIsomericSmiles(), structure_id:"ocl-" + await digest(nativeGraph),
      formula:formula.formula, molar_mass:formula.relativeWeight, architecture:"explicit_graph",
      net_charge:Array.from({length:count}, (_,i) => mol.getAtomCharge(i)).reduce((a,b)=>a+b,0)};
    let g = null, note = "Bu sonlu graf için konformer üretilemedi.";
    if (count <= 60) {
      try {
        const conformer = new O.ConformerGenerator(7301).getOneConformerAsMolecule(mol);
        if (conformer && conformer.getIDCode() === nativeGraph) {
          const atoms = Array.from({length:conformer.getAllAtoms()},(_,i)=>({
            element:conformer.getAtomLabel(i),x:conformer.getAtomX(i),y:conformer.getAtomY(i),z:conformer.getAtomZ(i)}));
          const bonds = Array.from({length:conformer.getAllBonds()},(_,i)=>({
            a:conformer.getBondAtom(0,i),b:conformer.getBondAtom(1,i),order:conformer.getBondOrder(i)}));
          g = usableGeometry({atoms,bonds,unit:"angstrom"});
          note = "OpenChemLib 9.26.0 ile üretilen tek bir izole başlangıç konformeri. Bağlar çizgiyle gösterilir. Sürükleyerek döndürebilirsin.";
        } else if (conformer) note = "Konformer oluşturulurken stereokimyasal kimlik değişti. Girdi korunuyor. Belirsiz merkezleri veya çift bağ geometrisini tanımlayarak yeniden dene.";
      } catch { /* The validated graph remains usable when 3D sampling fails. */ }
    } else note = "Tarayıcıda 3B üretim sınırı 60 açık atomdur. Bu grafın kaydı korunuyor. Daha büyük bir sonlu konformer için Yapı ve Konformer Enerjisi uygulamasını kullanabilirsin.";
    const raw = {schema:graphSchema,input_smiles:smiles,structure_id:structure.structure_id,engine:"OpenChemLib 9.26.0"};
    const exact = nativeGraph === O.Molecule.fromSmiles(known.structure.canonical_smiles).getIDCode();
    return {raw,structure,g,exact,note,native:false};
  }
  async function decode(data) {
    requireThat(data && typeof data === "object", "JSON bir yapı kaydı içermeli.");
    if (data.schema === nativeSchema || data.schema === identitySchema) return acceptNative(data);
    requireThat(data.schema === graphSchema, "Desteklenen dosyalar Yapı ve Konformer Enerjisi 1.0.0 kimlik ve sonuç kayıtları ile burada kaydedilen sonlu yapılardır.");
    const result = await buildGraph(data.input_smiles);
    requireThat(data.structure_id === result.structure.structure_id, "SMILES ile kaydedilmiş yapı kimliği eşleşmiyor.");
    return result;
  }
  function textLine(parent, tag, text) { const el = document.createElement(tag); el.textContent = text; parent.append(el); return el; }
  function label(key, value) { const dl = $("identity"); textLine(dl,"dt",key); textLine(dl,"dd",value); }
  function link(title, href) { const a=textLine($("next"),"a",title); a.href=href; }
  const kinds = {explicit_graph:"Girilen sonlu graf", linear_backbone:"Doğrusal ana zincir", star:"Yıldız"};
  const tactics = {as_drawn:"Girdide tanımlandığı gibi", syndiotactic:"Sindiyotaktik",isotactic:"İzotaktik",atactic:"Ataktik"};
  function show(result) {
    saved=result.raw; geo=result.g; yaw=.45; tilt=.32;
    const s=result.structure, dl=$("identity"), detail=$("detail"), next=$("next"); dl.replaceChildren();detail.replaceChildren();next.replaceChildren();
    label("Molekül formülü",s.formula);
    label("Modelin molar kütlesi",s.molar_mass.toLocaleString("tr-TR",{maximumFractionDigits:4})+" g/mol");
    label("Net yük",s.net_charge);label("Mimari",kinds[s.architecture] || s.architecture);
    if (Number.isInteger(s.dp)) label(s.mode==="star"?"Kol başına tekrar":"Modeldeki tekrar",s.dp);
    if (s.arm_count) label("Kol sayısı",s.arm_count);
    label("Uç grupları",s.h_caps == null ? "Sonlu SMILES içinde tanımlı" : s.h_caps+" H kapaması");
    if (s.tacticity) label("Taktisite",tactics[s.tacticity] || s.tacticity);
    if (s.stereo_pattern?.dyads) label("Doğrulanan ikililer",s.stereo_pattern.dyads);
    if (s.sequence) label("Birim sırası",s.sequence);
    if (result.native) label("Tanımsız stereo öğeleri",s.unspecified_stereo_count);
    textLine(detail,"code",s.canonical_smiles);
    textLine(detail,"p","Yapı kimliği"); textLine(detail,"code",s.structure_id);
    if (s.recipe_id) {textLine(detail,"p","Yapım reçetesi");textLine(detail,"code",JSON.stringify(s.recipe));textLine(detail,"code",s.recipe_id);}
    if (result.raw.artifact_id) {textLine(detail,"p","Kaynağın hesap kaydı");textLine(detail,"code",result.raw.artifact_id);}
    if (!result.native) {textLine(detail,"p","Girilen SMILES");textLine(detail,"code",result.raw.input_smiles);}
    $("geometry").hidden=!geo; $("geom-controls").hidden=!geo;
    $("geom-note").textContent=result.note || (geo?"İçe aktarılan ETKDGv3/MMFF94s taramasındaki en düşük enerjili kaydedilmiş izole konformer. Bağlar çizgiyle gösterilir. Sürükleyerek döndürebilirsin.":"Kaynak kayıtta 3B koordinat bulunmuyor.");
    $("scope").textContent = result.exact ?
      "Bu sonlu graf, su içindeki PEtOx trimer çiftinin model kimyasıyla eşleşiyor. Makroskopik faz hesabı için deneysel numunenin molar kütle dağılımı ve ölçüm koşulları ayrıca gerekli." :
      "Bu sonlu kimyasal yapı için eşleşmiş bir çözücü yörüngesi veya deneysel faz numunesi kayıtlı değil. Yapı bilgisi saklanabilir. Sonlu modelin kütlesi numunenin sayı veya ağırlık ortalamalı molar kütlesini belirlemez.";
    link("Bu yapı için çözeltideki zincir boyutunu incele", "chain.html?"+new URLSearchParams({source:"structure",system:"custom",name:s.formula+" · "+s.structure_id.slice(-10)}));
    if (result.exact) link("Bu kimyanın suyla temas filmini izle","contact.html");
    link("Bu yapı için çözelti viskozitesini incele", "viscosity.html?"+new URLSearchParams({source:"structure",system:"custom",name:s.formula+" · "+s.structure_id.slice(-10)}));
    const p=new URLSearchParams({system:"unsupported",name:(result.exact?"PEtOx trimeri · ":"")+s.formula+" · "+s.structure_id.slice(-10),structure_id:s.structure_id});
    link("Bu yapının faz verisi durumunu gör","phase.html?"+p);
    if (result.exact) link("PEtOx–su deneysel referans numunelerini aç","phase.html?system=petox");
    if (s.net_charge === 0 && Number.isFinite(s.molar_mass) && s.molar_mass > 0) {
      const condition = new URLSearchParams({finiteMass: String(s.molar_mass), formula: s.formula, charge: "0"});
      link("Bu yapı için seyreltik ozmotik sınırı incele", "ideal-osmotic.html?" + condition);
    }
    $("output").hidden=false; paint();
  }
  async function operation(task) {
    invalidate();const token=activeVersion;status("Yapı işleniyor…");await new Promise(resolve=>setTimeout(resolve,20));
    try {const result=await task(); if(token!==activeVersion)return;show(result);status("Yapı hazır.");}
    catch(error){if(token===activeVersion)status(error.message || "Kayıt açılamadı.",true);}
  }
  $("build").onclick=()=>operation(()=>buildGraph($("smiles").value.trim()));
  $("build-sketch").onclick=()=>operation(()=>{const g=sketch.getGraph();requireThat(g.finiteReady,"Uçları tanımlı, bağlı ve geçerli bir sonlu zincir çiz. Açık bağlantı noktaları için tekrar birimi uygulamasını kullan.");return buildGraph(g.smiles);});
  $("reference").onclick=()=>operation(()=>acceptNative(known));
  $("import-file").onchange=e=>{const file=e.target.files[0]; if(!file){invalidate();return;} operation(async()=>{requireThat(file.size<=10*1024*1024,"Dosya 10 MB sınırını aşıyor.");let obj;try {obj=JSON.parse(await file.text());}catch{throw new Error("Geçerli bir JSON dosyası seç.");}return decode(obj);});};
  $("save").onclick=()=>{if(!saved)return; const url=URL.createObjectURL(new Blob([JSON.stringify(saved,null,2)],{type:"application/json"}));const a=document.createElement("a");a.href=url;a.download="polimer-"+(saved.structure?.structure_id || saved.structure_id).slice(-12)+".json";document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),2000);};
  const colors={C:"#414d53",O:"#c34a41",N:"#345aa0",H:"#c7d0c5",F:"#417c5b",S:"#a17d24",Cl:"#3a805c",Br:"#9a5233",P:"#b66b35"};
  function paint() {
    if(!geo || $("output").hidden)return;
    const can=$("geometry"), w=can.clientWidth, h=can.clientHeight, dpr=window.devicePixelRatio || 1;
    if(!w||!h)return;
    can.width=Math.round(w*dpr);can.height=Math.round(h*dpr);
    const ctx=can.getContext("2d");ctx.scale(dpr,dpr);
    const indexes=geo.atoms.map((_,i)=>i).filter(i=>$("hs").checked||geo.atoms[i].element!=="H");
    if(!indexes.length)return;
    const center=["x","y","z"].map(key=>indexes.reduce((sum,i)=>sum+geo.atoms[i][key],0)/indexes.length);
    const r=Math.max(1,...indexes.map(i=>{const a=geo.atoms[i];return Math.hypot(a.x-center[0],a.y-center[1],a.z-center[2]);}));
    const scale=Math.min(w,h)*.43/r, ca=Math.cos(yaw),sa=Math.sin(yaw),ct=Math.cos(tilt),st=Math.sin(tilt);
    const pos=geo.atoms.map(a=>{const x=a.x-center[0],y=a.y-center[1],z=a.z-center[2],rz=-sa*x+ca*z;return [w/2+(ca*x+sa*z)*scale,h/2-(ct*y-st*rz)*scale,st*y+ct*rz];});
    const visible=new Set(indexes), elements=[];
    for(const b of geo.bonds)if(visible.has(b.a)&&visible.has(b.b))elements.push({depth:(pos[b.a][2]+pos[b.b][2])/2,b});
    for(const i of indexes)elements.push({depth:pos[i][2],i});
    elements.sort((a,b)=>a.depth-b.depth);ctx.lineCap="round";
    for(const elem of elements){if(elem.b){const a=pos[elem.b.a],b=pos[elem.b.b];ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.strokeStyle="#8b988e";ctx.lineWidth=Math.max(2,Math.min(5,scale*.11));ctx.stroke();}
      else{const p=pos[elem.i],e=geo.atoms[elem.i].element;ctx.beginPath();ctx.arc(p[0],p[1], Math.max(2.2,Math.min(13,scale*(e==="H"?.13:.23))),0,Math.PI*2);ctx.fillStyle=colors[e]||"#766092";ctx.fill();ctx.strokeStyle="#374c4799";ctx.lineWidth=.65;ctx.stroke();}}
  }
  let pointer=null;
  $("geometry").onpointerdown=e=>{pointer=[e.clientX,e.clientY];e.currentTarget.setPointerCapture(e.pointerId);};
  $("geometry").onpointermove=e=>{if(!pointer)return;yaw+=(e.clientX-pointer[0])*.014;tilt+=(e.clientY-pointer[1])*.014;pointer=[e.clientX,e.clientY];paint();};
  $("geometry").onpointerup=$("geometry").onpointercancel=()=>{pointer=null;};
  $("turn").onclick=()=>{yaw+=Math.PI/2;paint();};$("hs").onchange=paint;
  new ResizeObserver(paint).observe($("geometry"));
})();
