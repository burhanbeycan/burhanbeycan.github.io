(function () {
  "use strict";
  const byId = key => document.getElementById(key);
  const text = (id,value) => { byId(id).textContent=value; };
  const fmt = (value,digits=2) => Number(value).toLocaleString("tr-TR",{maximumFractionDigits:digits});
  const unit = value => value === null ? "—" : fmt(value) + " mPa·s";
  const model = window.HPAM_MODEL, engine=window.HPAM;
  if (!model || !engine) { text("loading","Model dosyası yüklenemedi."); return; }
  const svgNS="http://www.w3.org/2000/svg";
  function elSvg(tag, attrs={}, label) {
    const e = document.createElementNS(svgNS,tag);
    for (const [key,val] of Object.entries(attrs)) e.setAttribute(key,String(val));
    if (label !== undefined) e.textContent=label;
    return e;
  }
  function addRow(tbody, values) {
    const tr=document.createElement("tr");
    values.forEach(value => {const td=document.createElement("td");td.textContent=value;tr.append(td);});
    tbody.append(tr);
  }
  const conditions=model.test;
  const compositionList=[...new Set(conditions.map(d=>d.x.slice(0,2).join("|")))];
  const choose=byId("composition");
  choose.add(new Option("Serbest bileşim",""));
  compositionList.forEach(key=>{
    const [c,s]=key.split("|").map(Number);
    choose.add(new Option(fmt(c,0)+" ppm · "+fmt(s)+" g/L",key));
  });
  function summary(id, scores) {
    byId(id).replaceChildren();
    addRow(byId(id),["İkinci derece · en küçük kareler",fmt(scores.quadratic.mae,3),fmt(scores.quadratic.rmse,3)]);
    addRow(byId(id),["Kuvvet yasası",fmt(scores.baseline.mae,3),fmt(scores.baseline.rmse,3)]);
    addRow(byId(id),["İkinci dereceden Ridge",fmt(scores.ridge.mae,3),fmt(scores.ridge.rmse,3)]);
  }
  function draw(x, atComposition, range) {
    const svg=byId("curve"); svg.replaceChildren();
    const left=61,right=613,top=27,bottom=255;
    let data=[];
    if (range) for(let i=0;i<=90;i++){
      const g=range[0]+(range[1]-range[0])*i/90, point=[x[0],x[1],g];
      const b=engine.predict(model,point,"baseline"), r=engine.predict(model,point,"ridge"), q=engine.predict(model,point,"quadratic");
      if (r !== null && b !== null) data.push({g,r,b,q});
    }
    const yMax = Math.max(3, ...data.flatMap(d=>[d.r,d.b,d.q]),...atComposition.map(d=>d.y))*1.12;
    const toX = g=>left+(g-4)/52*(right-left);
    const toY = y=>bottom-y/yMax*(bottom-top);
    function line(x1,y1,x2,y2,color,more={}) {svg.append(elSvg("line",Object.assign({x1,y1,x2,y2,stroke:color},more)));}
    for(let i=0;i<=4;i++){
      const val=i*yMax/4, pos=toY(val);
      line(left,pos,right,pos,"#d6e0d9");
      svg.append(elSvg("text",{x:left-9,y:pos+4,"text-anchor":"end"},fmt(val,1)));
    }
    [4,10,20,30,40,50,56].forEach(g=>{
      const p=toX(g); line(p,bottom,p,bottom+5,"#587165");
      svg.append(elSvg("text",{x:p,y:bottom+20,"text-anchor":"middle"},fmt(g,0)));
    });
    svg.append(elSvg("text",{x:left,y:15},"η (mPa·s)"),elSvg("text",{x:(left+right)/2,y:304,"text-anchor":"middle"},"Kayma hızı (s⁻¹)"));
    if(!data.length) {
      svg.append(elSvg("text",{x:334,y:138,"text-anchor":"middle"},"Bu bileşim deney tasarımının dışında."));
      return;
    }
    for(const [key,color,dashed] of [["q","#644b78",false],["b","#98623b",true],["r","#19694f",true]]){
      const vals=data.map(d=>toX(d.g).toFixed(2)+","+toY(d[key]).toFixed(2)).join(" ");
      svg.append(elSvg("polyline",{points:vals,fill:"none",stroke:color,"stroke-width":2.8,"stroke-dasharray":dashed?"7 4":"none"}));
      if (Math.abs(data[0].g-data[data.length-1].g)<1e-7) svg.append(elSvg("circle",{cx:toX(data[0].g),cy:toY(data[0][key]),r:3,fill:color}));
    }
    if (engine.supported(model,x)) {
      line(toX(x[2]),top,toX(x[2]),bottom,"#7b8d81",{"stroke-dasharray":"3 4"});
      svg.append(elSvg("circle",{cx:toX(x[2]),cy:toY(engine.predict(model,x,"ridge")),r:5,fill:"#fff",stroke:"#19694f","stroke-width":2}));
    }
    atComposition.forEach(d=>svg.append(elSvg("circle",{cx:toX(d.x[2]),cy:toY(d.y),r:4,fill:"#253b30"})));
  }
  function render(){
    const x=["c","salt","shear"].map(id=>Number(byId(id).value));
    choose.value=compositionList.includes(x.slice(0,2).join("|"))?x.slice(0,2).join("|"):"";
    ["valC","valS","valG"].forEach((id,i)=>text(id,fmt(x[i])));
    const ok=engine.supported(model,x);
    const range=engine.shearInterval(model,x.slice(0,2));
    text("etaQ",unit(engine.predict(model,x,"quadratic")));text("etaB",unit(engine.predict(model,x,"baseline"))); text("etaR",unit(engine.predict(model,x,"ridge")));
    const subset=conditions.filter(r=>r.x[0]===x[0] && r.x[1]===x[1]);
    const obs=subset.find(r=>r.x[2]===x[2]);
    text("obs",obs?unit(obs.y):"—");
    text("rep",obs?(obs.n>1?obs.n+" tekrarın ortalaması · SS "+fmt(obs.sd)+" mPa·s":"Tek deney noktası"):"Bu koşulda ölçüm yayımlanmamış.");
    const msg=ok?(obs?"Seçili koşulda hem ölçüm hem tam veriyle eğitilmiş model gösteriliyor.":"Seçili koşul deney tasarımının içinde. En yakın ölçümler arasındaki model öngörüsü gösteriliyor."):"Seçili koşul deney tasarımının dışında. Derişimi, NaCl miktarını veya kayma hızını değiştirerek ölçülen bölgeye dönün.";
    text("domain",msg);byId("domain").className=ok?"reading rangeNote":"notice";
    text("supportText",range?"Bu bileşimde tasarımın kapsadığı kayma hızı "+fmt(range[0],1)+"–"+fmt(range[1],1)+" s⁻¹.":"Bu polimer ve tuz derişimlerinin birleşimi için ölçüm aralığı bulunmuyor.");
    draw(x,subset,range);
  }
  for(const id of ["c","salt","shear"])byId(id).addEventListener("input",render);
  choose.addEventListener("change",()=>{
    if (!choose.value) return;
    const [c,s]=choose.value.split("|").map(Number);
    const pts=conditions.filter(r=>r.x[0]===c && r.x[1]===s);
    const target=pts.reduce((a,b)=>Math.abs(a.x[2]-30)<=Math.abs(b.x[2]-30)?a:b);
    byId("c").value=c;byId("salt").value=s;byId("shear").value=target.x[2];render();
  });
  summary("summary",model.scores);
  text("validationInfo",model.nRuns+" yayımlanmış deney, "+conditions.length+" farklı koşul ve "+model.nGroups+" polimer–tuz bileşimi. Her turda bir bileşimin tüm kayma hızları ve tekrarları birlikte dışarıda tutuldu. Hata, her farklı koşulun ölçüm ortalamasından hesaplandı.");
  const names={baseline:"kuvvet yasası",quadratic:"ikinci dereceden en küçük kareler",ridge:"Ridge"};
  const best=Object.keys(names).sort((a,b)=>model.scores[a].rmse-model.scores[b].rmse)[0];
  text("conclusion","Bu deneyde en düşük RMSE değerini "+names[best]+" modeli verdi. Sonuç tek çalışmadaki dokuz HPAM–NaCl bileşimini kapsıyor.");
  for (const p of conditions) addRow(byId("validationRows"),[...p.x.map(v=>fmt(v)),p.n,fmt(p.y,3),fmt(p.baseline,3),fmt(p.quadratic,3),fmt(p.ridge,3)]);
  const foldChoice=byId("fold");
  model.folds.forEach((f,i)=>foldChoice.add(new Option(fmt(f.composition[0],0)+" ppm · "+fmt(f.composition[1])+" g/L",String(i))));
  function foldRender(){
    const f=model.folds[Number(foldChoice.value)];
    text("foldInfo",f.trainConditions+" eğitim ve "+f.testConditions+" sınama koşulu. Eğitim tasarım matrisinin rütbesi "+f.quadraticRank+"/10.");
    summary("foldScores",f.scores);
  }
  foldChoice.addEventListener("change",foldRender);foldRender();
  if (model.challenge.length) { byId("challengePanel").hidden=false;summary("challengeScores",model.challengeScores);text("challengeText",model.challenge.length+" koşul son modele uyumdan sonra sınandı.");}
  byId("save").addEventListener("click",()=>{
    const header=["source_runs","polymer_ppm","NaCl_g_L","shear_s-1","replicates","mean_eta_mPa_s","SD_mPa_s","OOF_power_mPa_s","OOF_quadratic_mPa_s","OOF_ridge_mPa_s"];
    const csv=[header,...conditions.map(r=>[r.ids.join("|"),...r.x,r.n,r.y,r.sd??"",r.baseline,r.quadratic,r.ridge])].map(row=>row.join(",")).join("\r\n");
    const href=URL.createObjectURL(new Blob([csv],{type:"text/csv;charset=utf-8"}));
    const a=document.createElement("a");a.href=href;a.download="hpam-formulasyon-sinamasi.csv";a.click();setTimeout(()=>URL.revokeObjectURL(href),1000);
  });
  byId("loading").hidden=true;byId("workspace").hidden=false;render();
})();
