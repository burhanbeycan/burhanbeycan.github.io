(function(){
"use strict";
const get=id=>document.getElementById(id);
const query=new URLSearchParams(location.search), NS="http://www.w3.org/2000/svg";
const sys=get("system"), sampleEl=get("sample"), fitEl=get("fit"), temp=get("temp"), phi=get("phi"), numPhi=get("phiValue");
const num=(v,d=3)=> Number(v).toLocaleString("tr-TR",{maximumFractionDigits:d});
const sci=(v,s=3)=>Math.abs(v)<1e-3 && v!==0 ? v.toExponential(s-1).replace(".",",").replace("e-"," ×10⁻").replace(/⁻(\d+)$/,(_,n)=>"⁻"+n.split("").map(d=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[d]).join("")) : num(v,Math.abs(v)<.1?6:4);
const theoretical={id:"constant",label:"Sabit integral χ",A:0,B:0,C:0};
let experiments=[], selectedSample=null, selectedFit=null, data=null, diagramCache=null;
function node(tag,text,className){const el=document.createElement(tag);if(text!=null)el.textContent=text;if(className)el.className=className;return el;}
function option(select,id,label){const o=node("option",label);o.value=id;select.add(o);}
function div(root,...children){root.replaceChildren(...children);}
function chooseSample(){selectedSample=data.samples.find(s=>s.id===sampleEl.value)||data.samples[0];experiments=data.measurements.filter(row=>(row.sampleId || row.series || row.sample) === selectedSample.id); }
function good(v){return typeof v==="number" && Number.isFinite(v);}
function boundsT(){
 if(sys.value==="theory") return [.45,.83];
 const ts=data.measurements.map(r=>r.T_K).filter(good);
 return ts.length?[Math.floor(Math.min(...ts)),Math.ceil(Math.max(...ts))]:[278,302];
}
function setSliderRange(keep){
 const [lo,hi]=boundsT();temp.min=lo;temp.max=hi;temp.step=sys.value==="theory"?.0001:.01;
 const target=sys.value==="theory"?.65:(experiments.length ? experiments.map(x=>x.T_K).sort((a,b)=>a-b)[Math.floor(experiments.length/2)] : (lo+hi)/2);
 temp.value=Math.max(lo,Math.min(hi,good(keep)?keep:target));
}
function initSystem(){
 const key=sys.value, isTheory=key==="theory", ps=window.EVIDENCE?.systems?.ps;
 const missingPS=key==="ps"&&!ps?.fit;
 const visible=key==="ps"&&!missingPS||isTheory;
 for (const k of ["sliders","phasePanel"]) get(k).hidden=!visible;
 get("sampleField").hidden=key!=="ps"||missingPS;
 get("fitField").hidden=key!=="ps"||missingPS;
 get("unavailable").hidden=visible||key==="petox";
 get("petoxPanel").hidden=key!=="petox";
 if(key==="petox"){if(typeof window.renderPetoxPanel==="function") window.renderPetoxPanel(get("petoxPanel"));return;}
 if(!visible){get("unknownText").textContent=key==="unsupported" ? "Aktarılan yapı: "+(query.get("name")||"Kullanıcı yapısı")+". Kimyasal graf ve mol kütlesi tek başına ölçülmüş bir etkileşim fonksiyonunu belirlemez." : "PS–siklohekzan deney tablosu bu oturumda yüklenemedi.";return;}
 if(key==="ps"){
  data=ps;
  sampleEl.replaceChildren();for(const s of ps.samples)option(sampleEl,s.id,s.label);
  const preferred=query.get("sample")||"ps-sip-123000";if(ps.samples.some(s=>s.id===preferred)) sampleEl.value=preferred;
  chooseSample();
  fitEl.replaceChildren();for(const [i,m] of (ps.fits||[ps.fit]).entries()) option(fitEl,i,m.label||m.id);
  const fits = ps.fits || [ps.fit];
  const requestedFit = fits.find(m => m.id === query.get("fit")) || ps.fit;
  fitEl.value = String(Math.max(0, fits.findIndex(m => m.id === requestedFit.id)));
  selectedFit=(ps.fits||[ps.fit])[+fitEl.value];
 }else{selectedSample={N:100};selectedFit=theoretical;experiments=[];}
 const pT=query.get("T");setSliderRange(key==="ps"&&pT&&Number.isFinite(+pT)?+pT:undefined);
 let p=query.get("phi");numPhi.value=phi.value=(p&&+p>0&&+p<1)?+p:.02;
 get("tLabel").textContent=isTheory?"Etkileşim parametresi χ":"Sıcaklık";
 get("diagramTitle").textContent=isTheory?"χ–φ faz diyagramı":"Sıcaklık–bileşim diyagramı";
 get("diagramSub").textContent=isTheory?"Sabit zincir hacim oranı N = 100. Eğriler teorik denge ve kararlılık sınırlarını gösterir.":"Daireler seçili numunenin ölçümleridir. Çizgiler seçilen termodinamik modelden hesaplanır.";
 get("obsLegend").hidden=isTheory;get("pointDetails").hidden=isTheory;diagramCache=null;update();
}
function chi(theta){return sys.value==="theory"?theta:selectedFit.A+selectedFit.B/theta;}
function rowPhi(r){return good(r.phi_p)?r.phi_p:undefined;}
function svgEl(tag,attrs={},text){const el=document.createElementNS(NS,tag);for(const[k,v]of Object.entries(attrs))el.setAttribute(k,v);if(text!==undefined)el.textContent=text;return el;}
function niceTicks(lo,hi,n){
 const raw=(hi-lo)/n;if(!(raw>0))return[lo];const pw=10**Math.floor(Math.log10(raw)),r=raw/pw;
 const step=pw*(r<=1?1:r<=2?2:r<=2.5?2.5:r<=5?5:10),t=[];
 for(let i=Math.ceil(lo/step-1e-10);i*step<=hi+step*1e-8&&t.length<30;i++) t.push(i*step);
 return t;
}
function tick(n){
 const a=Math.abs(n);if(a && (a<.001||a>50000))return n.toExponential(0).replace("e-","e−");
 return num(n,a<1?3:1);
}
function plane(svg,xmax,ylo,yhi,yTitle){
 const w=540,h=285,L=72,R=15,T=21,B=48,pw=w-L-R,ph=h-T-B;
 svg.setAttribute("viewBox",`0 0 ${w} ${h}`);svg.replaceChildren();
 const x=v=>L+v/xmax*pw,y=v=>T+(yhi-v)/(yhi-ylo)*ph;
 for(const v of niceTicks(0,xmax,5)){
   svg.append(svgEl("line",{x1:x(v),y1:T,x2:x(v),y2:h-B,stroke:"#e4e9e0"}),svgEl("text",{x:x(v),y:h-B+21,"text-anchor":"middle","font-size":14,fill:"#4c665a"},tick(v)));
 }
 for(const v of niceTicks(ylo,yhi,5)){
   svg.append(svgEl("line",{x1:L,y1:y(v),x2:w-R,y2:y(v),stroke:"#dbe3d8"}),svgEl("text",{x:L-8,y:y(v)+4,"text-anchor":"end","font-size":14,fill:"#4c665a"},tick(v)));
 }
 svg.append(svgEl("path",{d:`M${L} ${T}V${h-B}H${w-R}`,fill:"none",stroke:"#81988a"}),svgEl("text",{x:L+pw/2,y:h-4,"text-anchor":"middle","font-size":14,fill:"#294b3b"},"Polimer hacim kesri, φ"),svgEl("text",{x:14,y:T+ph/2,transform:`rotate(-90 14 ${T+ph/2})`,"text-anchor":"middle","font-size":14,fill:"#294b3b"},yTitle));
 const clipid="clip-"+svg.id,defs=svgEl("defs"),clip=svgEl("clipPath",{id:clipid});clip.append(svgEl("rect",{x:L,y:T,width:pw,height:ph}));defs.append(clip);svg.append(defs);
 const layer=svgEl("g",{"clip-path":"url(#"+clipid+")"});svg.append(layer);return{x,y,layer};
}
function path(g,pts,stroke,dash="",width=2.3){
 let d="",pen=false;for(const p of pts){
  if(!p || !p.every(Number.isFinite)){pen=false;continue;}
  d+=(pen?"L":"M")+g.x(p[0]).toFixed(2)+" "+g.y(p[1]).toFixed(2);pen=true;
 }
 g.layer.append(svgEl("path",{d,stroke,"stroke-width":width,fill:"none","stroke-dasharray":dash,"stroke-linejoin":"round"}));
}
function dot(g,px,py,color,r,title,outline){
 const d=svgEl("circle",{cx:g.x(px),cy:g.y(py),r,fill:color,stroke:outline||"#fff","stroke-width":1.6});if(title)d.append(svgEl("title",{},title));g.layer.append(d);return d;
}
function drawDiagram(p,T){
 const N=selectedSample.N,c1=selectedFit.C||0,[lo,hi]=boundsT();
 if(!diagramCache){
  const eqs=[];let ts=Array.from({length:67},(_,i)=>lo+(hi-lo)*i/66);
  const cc=FH.criticalChi(N,c1);const transition=sys.value==="theory"?cc:selectedFit.B/(cc-selectedFit.A);
  if(transition>lo&&transition<hi){ts.push(transition);ts.sort((a,b)=>a-b);}
  for(const t of ts){const e=FH.equilibrium(N,chi(t),c1),atCritical=Math.abs(t-transition)<1e-9;eqs.push([t,e,atCritical]);}
  diagramCache=eqs;
 }
 const xmax=Math.min(1,Math.max(.15,p*1.12,...diagramCache.map(([,e])=>(e.bin?.[1]||0)*1.03),...experiments.map(r=>(rowPhi(r)||0)*1.08)));
 const g=plane(get("diagram"),xmax,lo,hi,sys.value==="theory"?"Boyutsuz χ":"Sıcaklık (K)");
 for(const[k,color,dash]of [["bin","#166858",""],["spin","#af5d28","6 4"]])for(const j of [0,1]){
  path(g,diagramCache.map(([t,e,c])=>c?[e.pc,t]:e[k]?[e[k][j],t]:null),color,dash);
 }
 for(const [t,e,isCritical] of diagramCache)if(isCritical)dot(g,e.pc,t,"#242e29",4.3,"Model kritik noktası");
 for(const r of experiments){
   const x=rowPhi(r);if(!good(x))continue;
   const title="Ölçüm: φ = "+num(x,4)+" · "+num(r.T_K,2)+" K";
   const pt=dot(g,x,r.T_K,"#fff",3.7,title,"#415c95");pt.style.cursor="pointer";pt.setAttribute("aria-label",title);
   pt.onclick=()=>{temp.value=r.T_K;phi.value=x;numPhi.value=x;update();};
 }
 dot(g,p,T,"#b84d37",5.6,"Seçili bileşim");
}
function drawEnergy(p,eq){
 const [a,b]=eq.bin||[p,p];const z=eq.bin?b:p;
 const slope=eq.bin?eq.slope:FH.df(z,eq.N,eq.c0,eq.c1),base=FH.f(z,eq.N,eq.c0,eq.c1);
 const line=x=>base+(x-z)*slope,delta=x=>FH.f(x,eq.N,eq.c0,eq.c1)-line(x);
 const xMax=Math.min(1,Math.max(.035,p*1.18,b*1.13,eq.pc*1.25)),points=[];
 for(let i=0;i<=180;i++){const x=xMax*i/180;points.push([x,delta(x)]);}
 points.push([a,delta(a)],[b,delta(b)],[p,delta(p)]);points.sort((x,y)=>x[0]-y[0]);
 const scale=Math.max(1e-10,...points.map(q=>Math.abs(q[1]))),g=plane(get("energy"),xMax,-.075*scale,1.1*scale,"f − ℓ  (kBT / hücre)");
 path(g,points,"#166858","",2.9);
 const convex=points.map(([x,y])=>[x,eq.bin&&x>=a&&x<=b?0:y]);path(g,convex,"#845379","6 4",1.8);
 if(eq.bin){dot(g,a,0,"#845379",3.8,"Polimerce fakir faz");dot(g,b,0,"#845379",3.8,"Polimerce zengin faz");}
 dot(g,p,delta(p),"#b84d37",5.6,"Homojen karışım");
}
function basis(){
 const basis=get("basis"), sources=get("sources"),table=get("pointTable");
 if(sys.value==="theory"){
   div(basis,node("p","Bu ekran analitik bir yazılım kabul testidir. N = 100; χ sabittir. χ = 0,65 ve φ = 0,020 seçimi pozitif yerel eğriliğe rağmen iki denge fazı verir."),node("p","Beklenen binodal 0,009768 ve 0,271968; spinodal 0,038462 ve 0,200000.","muted"));
   sources.replaceChildren();return;
 }
 const s=selectedSample,m=selectedFit,note=node("p","Kafes hacmi oranı N = "+num(s.N,1)+". χ(φ,T) = "+num(m.A,5)+" "+(m.B>=0?"+ ":"− ")+num(Math.abs(m.B),3)+"/T"+(m.C?(m.C>=0?" + ":" − ")+num(Math.abs(m.C),5)+" φ":"")+". T kelvin cinsindendir.");
 const rows=experiments.filter(r=>good(r.T_K));
 div(basis,node("p",`${rows.length} yayımlanmış bulutlanma ölçümü bu numuneye ait. Bulutlanma bileşimleri dar dağılımlı ikili modelde yaklaşık binodal gözlemleri olarak kullanılır. Spinodal çizgisi modelden türetilir.`),note);
 if(s.dispersity)basis.append(node("p","Molar kütle ve dağılım: Mw = "+num(s.Mw_g_mol,0)+" g/mol, Đ = "+num(s.dispersity,3)+".","muted"));
 const split=s.split==="holdout"?"Yayın bazında dışarıda bırakılan sınama numunesi.":"Kalibrasyonda kullanılan numune.";
 basis.append(node("p",split+" Numuneye ait ölçüm aralığı "+num(s.T_range_K[0],2)+"–"+num(s.T_range_K[1],2)+" K ve φ = "+num(s.phi_range[0],5)+"–"+num(s.phi_range[1],5)+". 0,1 MPa.","muted"));
 basis.append(node("p","N = (Mw/Đ) v̄PS / Vs. Referans değerler v̄PS = 0,934 mL/g ve Vs = 108,04 mL/mol. Dar dağılımlı zincirler tek bir polimer türüyle temsil edilir.","muted"));
 const criticalIntercept=FH.criticalChi(s.N,m.C||0),Tc=m.B/(criticalIntercept-m.A);
 basis.append(node("p","Bu numunenin model kritik noktası "+num(Tc,2)+" K, φc = "+num(FH.criticalPhi(s.N,m.C||0),5)+"."));
 if(m.evalMetric)basis.append(node("p","Ayrı yayının "+m.evalMetric.n+" bulutlanma kaydında RMSE = "+num(m.evalMetric.rmse_K,2)+" K; kalibrasyonda "+m.trainMetric.n+" kayıt. RMSE bu kayıtlardaki ortalama hata büyüklüğünü özetler."));
 const check=node("a","Kalibrasyon ve dış sınama sonuçlarını incele");check.href="validation.html";basis.append(check);
 if(m.validationText)basis.append(node("p",m.validationText));
 const ul=node("ul",undefined,"sourceList");
 for(const ref of data.sources||[]){const url=ref.url||(ref.doi?"https://doi.org/"+ref.doi:null);if(!url)continue;const l=node("li"),a=node("a",ref.name||ref.title||ref.doi||"Ölçüm kaynağı");a.href=url;a.target="_blank";a.rel="noopener";l.append(a);ul.append(l);}
 const csv=node("a","97 ölçüm ve iki modelin tahminleri · CSV");csv.href="ps-cloud-points.csv";csv.download="ps-cloud-points.csv";div(sources,ul,csv);
 const wrap=node("div",undefined,"tableWrap"),tbl=node("table"),head=tbl.createTHead().insertRow();for(const label of ["T (K)","φ","Ölçüm","Seç"])head.append(node("th",label));
 const body=tbl.createTBody();for(const r of rows){const tr=body.insertRow();tr.insertCell().textContent=num(r.T_K,2);const x=rowPhi(r);tr.insertCell().textContent=good(x)?num(x,5):"—";tr.insertCell().textContent="Bulutlanma";const cell=tr.insertCell();if(good(x)){const b=node("button","İncele");b.onclick=()=>{temp.value=r.T_K;phi.value=x;numPhi.value=x;update();};cell.append(b);}}
 wrap.append(tbl);div(table,wrap);
}
function update(){
 if(get("phasePanel").hidden)return;
 const p=Math.max(.00001,Math.min(.99999,Number(numPhi.value))),T=+temp.value,N=selectedSample.N,c0=chi(T),C=selectedFit.C||0,eq=FH.equilibrium(N,c0,C),res=FH.classify(p,eq);
 const names={one:"Küçük dalgalanma söner",meta:"Küçük dalgalanma söner",unstable:"Küçük dalgalanma büyür",critical:"Kritik bileşim",spinodal:"Spinodal sınırı",boundary:"Binodal sınırı"};
 const details={one:"Kararlı tek faz",meta:"Metastabil homojen karışım",unstable:"Spinodal kararsızlık",critical:"Faz bileşimleri çakışır",spinodal:"Yerel eğrilik sıfır",boundary:"Diğer fazın payı sıfır"};
 get("phaseCount").textContent=res.count===2?"İki faz":"Tek faz";get("phaseDetail").textContent=details[res.label];get("local").textContent=names[res.label];get("curvature").textContent="f″ = "+(res.curvature>=0?"+":"")+sci(res.curvature);
 get("lever").textContent=res.richFraction==null?"—":"%"+num(100*res.richFraction,2);
 get("leverNote").textContent=res.richFraction==null?"Eşfaz oranı uygulanmaz":"Toplam karışım hacmindeki pay";
 get("tUnit").textContent=sys.value==="theory"?num(T,4):num(T,2)+" K";
 const rich=res.count===2?" Faz bileşimleri φ₁ "+(eq.diluteUnderflow?"< 5 ×10⁻³²⁴":"= "+sci(eq.bin[0]))+" ve φ₂ = "+sci(eq.bin[1])+".":"";
 const note={one:"Homojen bileşim bu modelde küresel serbest enerji minimumundadır.",meta:"İki faz daha düşük serbest enerji verir. Homojen durum küçük dalgalanmalara yerel olarak dayanır; ayrışma için yeterli bir çekirdek gerekir.",unstable:"Homojen durum negatif eğriliktedir. Yeterince uzun dalga boylu küçük dalgalanmalar büyüyebilir.",critical:"Modelin kritik noktası yerel ikinci ve üçüncü türevlerin sıfıra gittiği bileşimdir.",spinodal:"Seçili bileşim homojen kararlılık sınırında yer alır.",boundary:"Seçili bileşim denge fazlarından birinin bileşimine eşittir."};const beyond=sys.value==="ps" && (p<selectedSample.phi_range[0]||p>selectedSample.phi_range[1]||T<selectedSample.T_range_K[0]||T>selectedSample.T_range_K[1]);
 get("stateNote").textContent=note[res.label]+rich+(beyond?" Seçim, bu numunede ölçülen sıcaklık veya bileşim aralığının dışında.":"");
 drawDiagram(p,T);drawEnergy(p,eq);
 const queryK=new URLSearchParams({N:String(N),chi:String(c0),chi1:String(C),phi:String(p),mode:"quench"});if(sys.value==="ps")queryK.set("T",T);
 const osmoticLink = get("osmoticTransfer");
 osmoticLink.hidden = sys.value !== "ps";
 if (sys.value === "ps") {
   const condition = new URLSearchParams({sample: selectedSample.id, T: String(T), phi: String(p), model: selectedFit.id === "fh0" ? "fh0" : "cloud"});
   osmoticLink.href = "osmotic.html?" + condition;
 }
 const xfer=get("transfer");xfer.hidden=false;xfer.href="kinetics.html?"+queryK;
 get("transferText").textContent="Zincir hacim oranı, χ fonksiyonu ve ortalama bileşim aynı yerel serbest enerjiyle faz alanı çözücüsüne aktarılır. Kinetik zaman ve uzunluk indirgenmiş birimlerdedir.";
}
sys.addEventListener("change",()=>{initSystem();if(!get("phasePanel").hidden)basis();});
sampleEl.addEventListener("change",()=>{chooseSample();diagramCache=null;update();basis();});
fitEl.addEventListener("change",()=>{selectedFit=(data.fits||[data.fit])[+fitEl.value];diagramCache=null;update();basis();});
temp.addEventListener("input",update);phi.addEventListener("input",()=>{numPhi.value=phi.value;update();});
numPhi.addEventListener("input",()=>{if(numPhi.value!==""&&numPhi.validity.valid&&+numPhi.value>0&&+numPhi.value<1){phi.value=numPhi.value;update();}});numPhi.addEventListener("change",()=>{if(!numPhi.validity.valid||numPhi.value===""){numPhi.value=phi.value;update();}});
const requested=query.get("system");if(requested==="unsupported")option(sys,"unsupported","Aktarılan özel yapı");if(["ps","petox","theory","unsupported"].includes(requested))sys.value=requested;
initSystem();if(!get("phasePanel").hidden)basis();
})();
