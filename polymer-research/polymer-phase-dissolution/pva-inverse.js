'use strict';
const readout = document.getElementById.bind(document);
const research = globalThis.PvaResearch;
const fmt = (x, digits=1) => Number(x).toLocaleString('tr-TR', {maximumFractionDigits:digits});
const flow = (c,t,g) => globalThis.PvaInverseModel.predict(research.flow,c,t,g);
const inputs = id => readout(id).value === '' ? NaN : Number(readout(id).value);
function interpolateCurve(curve, rate) {
  if (!curve || !(rate>0)) return null;
  const exact=curve.find(p=>Math.abs(Math.log(p[0]/rate))<1e-8);
  if(exact) return exact[1];
  const hi=curve.findIndex(p=>p[0]>rate);
  if(hi<=0) return null;
  const [x0,y0]=curve[hi-1], [x1,y1]=curve[hi];
  if(Math.log10(x1/x0)>0.35) return null;
  return y0*Math.exp(Math.log(y1/y0)*Math.log(rate/x0)/Math.log(x1/x0));
}
function designMass(mn, dh, c, t, rate) {
  const all=Object.values(research.mass.grades).filter(s=>s.dh_mol_pct===dh).sort((a,b)=>a.mn_g_mol-b.mn_g_mol);
  const a=all.filter(s=>s.mn_g_mol<=mn).pop(), b=all.find(s=>s.mn_g_mol>=mn);
  const miss = reason => ({value:null,reason});
  if (!a || !b) return miss('Seçilen hidroliz derecesinde bu Mₙ değerini çevreleyen iki doğrulanmış numune yok.');
  const k=Number(c).toFixed(1);
  const va=interpolateCurve(a.curves[k],rate), vb=interpolateCurve(b.curves[k],rate);
  if(va===null || vb===null) return miss('Bu derişim ve hızda gerekli numunelerin yakın ölçüm noktaları yok. Kütlece %4 ve 60 s⁻¹ koşulunu inceleyebilirsiniz.');
  let eta=va;
  if(a!==b) eta=va*Math.exp(Math.log(vb/va)*Math.log(mn/a.mn_g_mol)/Math.log(b.mn_g_mol/a.mn_g_mol));
  if(t===298) {
    const target=flow(c,298,rate), ref=flow(c,298.15,rate);
    if(!target.supported || !ref.supported) return miss('Bu koşul için 298,00 K sıcaklık düzeltmesi desteklenmiyor. 298,15 K deney sıcaklığı seçilebilir.');
    eta*=target.eta_cP/ref.eta_cP;
  }
  return {value:eta,reason:(a===b?'Tek numunenin sayısallaştırılmış akış eğrisi.':'Mₙ referansları '+fmt(a.mn_g_mol,0)+' ve '+fmt(b.mn_g_mol,0)+' g/mol. Đ değerleri '+fmt(a.dispersity,2)+' ve '+fmt(b.dispersity,2)+'.')+(t===298?' Sıcaklık oranı başka bir PVA ürününden aktarıldı.':' Deney sıcaklığında hesaplandı.')};
}
const xgrid=(a,b,n)=>Array.from({length:n},(_,i)=>a*(b/a)**(i/(n-1)));
function result() {
  const mode=readout('mode').value;
  if(mode==='transfer'){
    const idx=Number(readout('curveType').value);
    const curve=research.transfer[idx], selected=Number(readout('transferRate').value);
    const p=curve.rows[selected]||curve.rows[0];
    return {value:p[2],kind:'polyBERT + kimya · aile dışı tahmin',conditions:curve.material+' · '+fmt(curve.t,2)+' K\n'+fmt(p[0])+' s⁻¹',
      reason:'Yayımlanmış viskozite '+fmt(p[1],2)+' cP. Bu noktada mutlak yüzde hata %'+fmt(100*Math.abs(p[2]/p[1]-1))+'.',
      legend:'Yeşil · aile eğitimden çıkarıldığında hesaplanan tahminler. Turuncu · gerçek ölçümler.',
      descr:'Tekrar birimi gösterimi ve derişim birimi protokole özgüdür. Aileler arasında numune ve laboratuvar da değişebilir. Kaynak: '+curve.doi+'.',
      x:p[0],axis:'Kayma hızı (s⁻¹)',points:curve.rows.map(q=>[q[0],q[2]]),measured:curve.rows.map(q=>[q[0],q[1]])};
  }
  if(mode==='flow') {
    const c=inputs('flowC'),t=inputs('flowT'),g=inputs('flowRate');
    const r=flow(c,t,g);
    return {value:r.eta_cP,kind:'Görünür viskozite · ters model',conditions:fmt(t,2)+' K · %'+fmt(c)+' PVA\n'+fmt(g)+' s⁻¹',reason:r.supported?'Merck 141351 için parametreleri akış eğrilerinden öğrenilen model. Bu ürünün Mₙ değeri kaynakta raporlanmadı.':'Seçili koşul ortak ölçüm alanının dışında. Desteklenen komşu eğriler gerekli.',legend:'Yeşil çizgi · ters model. Nokta · seçilen hız.',descr:'Çözücü: deiyonize su. Veri kaynağı: 16 derişim–sıcaklık eğrisi, 244 reometre ölçümü.',x:g,axis:'Kayma hızı (s⁻¹)',points:xgrid(.12,155,145).map(x=>[x,flow(c,t,x).eta_cP])};
  }
  if(mode==='mn') {
    const m=inputs('mn'),dh=inputs('dh'),c=inputs('gradeC'),t=inputs('gradeTemp'),g=inputs('gradeRate');
    const r=designMass(m,dh,c,t,g);
    const h=dh/100,M0=h*44.053+(1-h)*86.090;
    return {value:r.value,kind:'Teorik Mₙ · araştırma tahmini',conditions:fmt(t,2)+' K · %'+fmt(c)+' PVA\n'+fmt(g)+' s⁻¹',reason:r.reason,legend:'Yeşil çizgi · seçili teorik Mₙ için desteklenen hızlar. Nokta · seçilen hız.',descr:'Ortalama tekrar birimi '+fmt(M0,2)+' g/mol · kimyasal OH yoğunluğu '+fmt(1000*h/M0,2)+' mmol/g. Bunlar hidroliz derecesinden hesaplanan kuru zincir değerleridir.',x:g,axis:'Kayma hızı (s⁻¹)',points:xgrid(.006,1200,350).map(x=>[x,designMass(m,dh,c,t,x).value])};
  }
  const id=readout('material').value, ages=research.aging.time_predictions.filter(p=>p.material===id);
  const s=ages[0], hours=inputs('hours');
  const estimate=h=>s.fresh_cP*(1+research.aging.K_held_two_materials*s.wt_percent**2*h);
  const ok=Number.isFinite(hours)&&hours>=0&&hours<=16;
  return {value:ok?estimate(hours):null,kind:'Düşen bilye · koşullu yaşlanma tahmini',conditions:'303,15 K · %'+fmt(s.wt_percent,2)+' PVA\n'+fmt(hours)+' saat',reason:ok?'Başlangıç viskozitesi bu numunenin ölçümünden gelir. Erken zaman yasası diğer numunelerin yaşlanma hızlarıyla kalibre edildi.':'Bu çalışma için erken zaman alanı 0–16 saattir.',legend:'Yeşil çizgi · 0–16 saat kinetik hesabı. Turuncu noktalar · sayısallaştırılmış ölçümler.',descr:'Numune '+id+' · başlangıç '+fmt(s.fresh_cP)+' cP · hız sabiti '+fmt(research.aging.K_held_two_materials*s.wt_percent**2,4)+' saat⁻¹. Düşen bilye ölçümleri için tanımlı bir sürekli kayma hızı raporlanmadı.',x:hours,axis:'Bekleme süresi (saat)',linear:true,points:Array.from({length:33},(_,i)=>[i/2,estimate(i/2)]),measured:ages.filter(p=>p.age_h<=16).map(p=>[p.age_h,p.observed_cP]).concat([[0,s.fresh_cP]])};
}
function graph(r) {
  const canvas=readout('plot'),box=canvas.getBoundingClientRect(),w=box.width,h=box.height;
  if(!w||!h)return;
  const dpr=Math.min(devicePixelRatio||1,3);canvas.width=w*dpr;canvas.height=h*dpr;
  const ctx=canvas.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
  const ar=r.points.filter(p=>Number.isFinite(p[1])&&p[1]>0);
  if(!ar.length){ctx.fillStyle='#486366';ctx.font='14px system-ui';ctx.fillText('Desteklenen eğri yok',22,h/2);return;}
  const pts=ar.concat(r.measured||[]),fx=x=>r.linear?x:Math.log10(x),fy=y=>r.linear?y:Math.log10(y);
  const x0=Math.min(...pts.map(p=>fx(p[0]))),x1=Math.max(...pts.map(p=>fx(p[0])));
  let y0=Math.min(...pts.map(p=>fy(p[1]))),y1=Math.max(...pts.map(p=>fy(p[1])));const dy=(y1-y0||1)*.1;y0-=dy;y1+=dy;
  const left=66,right=w-12,top=14,bottom=h-43,X=x=>x0===x1?(left+right)/2:left+(fx(x)-x0)/(x1-x0)*(right-left),Y=y=>bottom-(fy(y)-y0)/(y1-y0)*(bottom-top);
  ctx.font='12px system-ui';ctx.fillStyle='#466062';ctx.lineWidth=1;
  for(let i=0;i<=3;i++){const z=y0+(y1-y0)*i/3,y=r.linear?z:10**z;ctx.strokeStyle='#dbe4df';ctx.beginPath();ctx.moveTo(left,Y(y));ctx.lineTo(right,Y(y));ctx.stroke();ctx.textAlign='right';ctx.fillText(fmt(y,y<10?2:0),left-6,Y(y)+3);}
  for(let i=0;i<=(x0===x1?0:3);i++){const z=x0+(x1-x0)*i/3,x=r.linear?z:10**z;ctx.textAlign=x0===x1?'center':i===0?'left':i===3?'right':'center';ctx.fillText(fmt(x,x<10?2:0),X(x),bottom+16);}
  ctx.textAlign='center';ctx.fillText(r.axis,(left+right)/2,h-4);ctx.textAlign='left';ctx.fillText('η (cP)',3,11);
  let prev=false;ctx.beginPath();ctx.strokeStyle='#08735d';ctx.lineWidth=2.5;
  for(const [x,y] of r.points){if(!(y>0)){prev=false;continue;}if(prev)ctx.lineTo(X(x),Y(y));else ctx.moveTo(X(x),Y(y));prev=true;}ctx.stroke();
  for(const [x,y] of r.measured||[]){ctx.fillStyle='#ab5729';ctx.beginPath();ctx.arc(X(x),Y(y),3.5,0,2*Math.PI);ctx.fill();}
  if(r.value>0&&Number.isFinite(r.x)&&fx(r.x)>=x0&&fx(r.x)<=x1){ctx.strokeStyle='#113e36';ctx.fillStyle='white';ctx.lineWidth=2;ctx.beginPath();ctx.arc(X(r.x),Y(r.value),5,0,2*Math.PI);ctx.fill();ctx.stroke();}
}
function update() {
 const mode=readout('mode').value;
 for(const block of document.querySelectorAll('[data-mode]'))block.hidden=block.dataset.mode!==mode;
 const r=result();
 for(const [id,text] of Object.entries({kind:r.kind,value:r.value>0?fmt(r.value,r.value<10?2:1)+' cP':'—',message:r.reason,conditions:r.conditions,legend:r.legend,descriptor:r.descr}))readout(id).textContent=text;
 graph(r);
}
for(const id of [...new Set(research.aging.time_predictions.map(r=>r.material))]){const op=document.createElement('option');op.value=id;op.textContent=id;readout('material').appendChild(op);}
const names={PVA:'PVA',HPAM:'Hidrolize poliakrilamid',CMC:'Sodyum karboksimetilselüloz',ALG:'Sodyum aljinat'};
const option=(sel,v,txt)=>{const op=document.createElement('option');op.value=v;op.textContent=txt;sel.appendChild(op);};
for(const m of [...new Set(research.transfer.map(c=>c.material))]) option(readout('chemType'),m,names[m]||m);
function syncRates(){
 const s=readout('transferRate');s.replaceChildren();
 const curve=research.transfer[Number(readout('curveType').value)];curve.rows.forEach((r,i)=>option(s,i,fmt(r[0],2)));
}
function syncCurves(){
 const s=readout('curveType');s.replaceChildren();
 research.transfer.forEach((r,i)=>{if(r.material===readout('chemType').value)option(s,i,r.label);});syncRates();
}
readout('chemType').addEventListener('input',syncCurves);readout('curveType').addEventListener('input',syncRates);syncCurves();
for(const row of research.evidence){const tr=document.createElement('tr');for(const value of row){const td=document.createElement('td');td.textContent=value;tr.appendChild(td);}readout('evidence').appendChild(tr);}
readout('inputs').addEventListener('input',update);readout('mode').addEventListener('change',update);window.addEventListener('resize',update);update();
