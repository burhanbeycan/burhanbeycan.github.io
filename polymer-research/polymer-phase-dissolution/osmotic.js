(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  const sys = window.EVIDENCE?.systems?.ps;
  const report = window.TRANSFER_RESULTS;
  const external = window.OSMOTIC_PRESSURE_DATA;
  const query = new URLSearchParams(location.search);
  const format = (v, significant = 4) => v.toLocaleString('tr-TR', {maximumSignificantDigits: significant});
  const numberInput = value => String(Number(value.toPrecision(12)));
  function scientific(x) {
    if (x === 0) return '0';
    if (Math.abs(x) >= .001 && Math.abs(x) < 100000) return format(x);
    const [coefficient, exponent] = x.toExponential(2).split('e');
    const superscripts = {'-':'⁻','+':'','0':'⁰','1':'¹','2':'²','3':'³','4':'⁴','5':'⁵','6':'⁶','7':'⁷','8':'⁸','9':'⁹'};
    return format(Number(coefficient)) + ' × 10' + [...exponent].map(k => superscripts[k]).join('');
  }
  function halt(message) {
    $('controls').hidden = true;
    $('results').hidden = true; $('external').hidden=true;
    $('unavailable').hidden = false;
    $('unavailable').textContent = message + ' ';
    const a=document.createElement('a'); a.href='ideal-osmotic.html'; a.textContent='Kendi numunemin seyreltik sınırını incele'; $('unavailable').append(a);
  }
  if (!sys?.samples?.length || !sys?.fits || !window.FH || !window.Osmosis) {
    halt('Sayısal model veya numune verileri yüklenemedi.');
    return;
  }
  const allSamples = external ? [external,...sys.samples] : sys.samples;
  const sampleId = query.get('sample');
  if ((query.has('system') && query.get('system') !== 'ps') || (sampleId && !allSamples.some(s => s.id === sampleId))) {
    halt('Aktarılan kimya veya numune için doğrulanmış bir PS–siklohekzan serisi bulunamadı.');
    return;
  }
  const {polymer_specific_volume_mL_g:v, solvent_molar_volume_mL_mol:Vs} = sys.volume_model;
  const models = new Map();
  const addModel = (key,label,model) => {
    if (!model || ![model.A,model.B,model.C].every(Number.isFinite)) return;
    models.set(key, model); $('model').add(new Option(label,key));
  };
  addModel('joint','Bulutlanma + ozmotik A₂',report?.models?.find(m => m.id === 'combined'));
  addModel('cloud','Bulutlanma · χ(φ,T)',report?.models?.find(m => m.id === 'cloud') || sys.fits.find(m => m.id === 'fh1'));
  addModel('fh0','Bulutlanma · χ(T)',sys.fits.find(m => m.id === 'fh0'));
  for (const s of allSamples) $('sample').add(new Option(s.label,s.id));
  $('sample').value = sampleId || (external ? external.id : 'ps-sip-123000');
  if (models.has(query.get('model'))) $('model').value = query.get('model');
  $('temperature').value = query.has('T') ? query.get('T') : (external && !sampleId ? '307.15' : '298.15');
  const cMax = .99999 * 1000/v;
  $('concentration').max = String(cMax);
  $('concentration').value = query.has('c') ? query.get('c') : query.has('phi') ? numberInput(Number(query.get('phi'))*1000/v) : (external && !sampleId ? '20.1' : '5');
  const ns='http://www.w3.org/2000/svg';
  function svg(tag,attributes,label) {
    const item = document.createElementNS(ns,tag);
    Object.entries(attributes).forEach(([k,val]) => item.setAttribute(k,val));
    if(label != null) item.textContent=label;
    return item;
  }
  const inputIds=['sample','model','temperature','concentration'];
  let latest = null, chartRows = [], scanLimit=0;
  function invalid(message) {
    $('inputError').textContent=message; $('inputError').hidden=false; $('results').hidden=true;latest=null;
  }
  function makeGrid(s,selected,max) {
    const cs=Array.from({length:121},(_,i) => max * i / 120);
    for (const p of s.eq.bin || []) {
      const c=p*1000/v; if (c>0 && c<max) cs.push(c);
    }
    cs.push(selected);
    return [...new Set(cs)].sort((a,b)=>a-b).map(c => Osmosis.evaluate(c,s));
  }
  function plot(s,point,max) {
    const plot=$('plot'), isPi=$('axis').value==='pressure';
    chartRows=makeGrid(s,point.c,max); plot.replaceChildren();
    const left=88,right=618,top=15,bottom=234;
    const ideal = r => isPi ? r.ideal : Osmosis.R*s.T/s.Mn;
    const measure = r => isPi ? r.pressure : r.reduced;
    const showObservations=external && $('sample').value===external.id ? external.observations.filter(r=>Math.abs(r.T_K-s.T)<1e-8 && r.c_g_L<=max) : [];
    const observedValue=r=>isPi?r.Pi_kPa:r.Pi_kPa/r.c_g_L;
    $('pointKey').hidden=!showObservations.length;
    const ys=[...chartRows.flatMap(r => [measure(r),ideal(r)]), ...showObservations.map(observedValue)];
    let min=isPi?0:Math.min(...ys),upper=Math.max(...ys);
    const range=Math.max(upper-min,Math.abs(upper)*.08,1e-10);
    min=Math.max(0,min-range*.08);upper+=range*.09;
    const X=c=>left+(right-left)*c/max;
    const Y=p=>bottom-(bottom-top)*(p-min)/(upper-min);
    const bin=s.eq.bin, low=bin?bin[0]*1000/v:0, high=bin?bin[1]*1000/v:0;
    const hasBand=!!bin && low<max;
    $('band').hidden=!hasBand;
    if (hasBand) plot.append(svg('rect',{x:X(low),y:top,width:X(Math.min(max,high))-X(low),height:bottom-top,fill:'#e4efdf'}));
    const line=(x1,y1,x2,y2,color='#c9d6cb')=>plot.append(svg('line',{x1,y1,x2,y2,stroke:color,'stroke-width':1}));
    const text=(x,y,msg,anchor='middle',fontSize=15)=>plot.append(svg('text',{x,y,'font-size':fontSize,fill:'#3b5749','text-anchor':anchor},msg));
    for(let i=0;i<=4;i++){const value=min+(upper-min)*i/4;line(left,Y(value),right,Y(value));text(left-8,Y(value)+5,scientific(Number(value.toPrecision(3))),'end',14);}
    line(left,top,left,bottom,'#83998b');line(left,bottom,right,bottom,'#83998b');
    for(let i=0;i<=4;i++){const val=max*i/4;line(X(val),bottom,X(val),bottom+4,'#83998b');text(X(val),bottom+22,format(val,3));}
    text((left+right)/2,285,'Polimer derişimi c (g/L)');
    const yTitle=svg('text',{x:0,y:0,transform:'translate(18 '+(top+bottom)/2+') rotate(-90)','text-anchor':'middle','font-size':14,fill:'#284c3f'},isPi?'Π (kPa)':'Π / c (kPa L/g)');plot.append(yTitle);
    function trace(get,color,dash) {
      const d=chartRows.map((r,i)=>(i?'L':'M')+X(r.c).toFixed(4)+','+Y(get(r)).toFixed(4)).join(' ');
      const at={d,fill:'none',stroke:color,'stroke-width':2.7,'stroke-linejoin':'round'};if(dash)at['stroke-dasharray']='6 4';plot.append(svg('path',at));
    }
    trace(ideal,'#8c664b',true);trace(measure,'#176c54',false);
    for(const row of showObservations){
      const observed=svg('circle',{cx:X(row.c_g_L),cy:Y(observedValue(row)),r:4.6,fill:'white',stroke:'#775795','stroke-width':2.1});
      observed.append(svg('title',{},format(row.c_g_L)+' g/L · '+format(row.Pi_kPa)+' kPa'));
      plot.append(observed);
    }
    plot.append(svg('circle',{cx:X(point.c),cy:Y(measure(point)),r:5,fill:'#141f23',stroke:'white','stroke-width':1.4}));
    plot.setAttribute('aria-label',isPi?'Çözelti derişimine göre modelden hesaplanan denge ozmotik basıncı':'Çözelti derişimine göre denge basıncının polimer derişimine oranı');
    $('figureNote').textContent=isPi?'Yeşil eğri, iki fazlı bölgede ortak teğetten hesaplanan denge basıncını izler. İşaretli nokta seçtiğiniz derişimdir.':'Sonsuz seyreltmede Π/c, RT/Mₙ değerine yaklaşır. İki fazlı bölgede sabit Π, artan toplam c’ye bölündüğü için Π/c azalır.';
  }
  function update(origin) {
    const T=$('temperature').valueAsNumber,c=$('concentration').valueAsNumber;
    if(!Number.isFinite(T)||!$('temperature').validity.valid||!Number.isFinite(c)||c<0||c>cMax){invalid('278–318,15 K arasında sıcaklık ve 0–'+format(cMax,6)+' g/L arasında derişim girin.');return;}
    $('inputError').hidden=true; $('results').hidden=false;
    const sample=allSamples.find(x=>x.id===$('sample').value),fit=models.get($('model').value);
    const s=Osmosis.state(T,sample.Mn_g_mol,v,Vs,fit),point=Osmosis.evaluate(c,s);
    latest={s,point,sample,model:$('model').value};
    if(origin!=='scan') {
      const target=c<=15?30:c<=50?100:c<=125?250:c<=250?500:cMax;
      const obsCap=external && sample.id===external.id && external.observations.some(row=>Math.abs(row.T_K-T)<1e-8)?270:0;
      scanLimit=Math.max(c,target,obsCap);
      $('scan').max=String(scanLimit);
    }
    $('scan').value=String(c);
    $('pi').textContent=scientific(point.pressure)+' kPa';
    $('ideal').textContent=scientific(point.ideal)+' kPa';
    $('a2').textContent=scientific(s.A2);
    $('phase').textContent=point.phase.count===2?'İki faz dengesi · basınç platosu':(point.phase.label==='boundary'?'Denge fazının sınırı':point.phase.label==='critical'?'Modelin kritik noktası':'Tek fazlı denge');
    $('mn').textContent=format(sample.Mn_g_mol/1000,5)+' kg/mol';$('disp').textContent=sample.dispersity==null?'Bildirilmemiş':format(sample.dispersity);$('phi').textContent=scientific(point.phi);
    $('plotCondition').textContent=format(T,6)+' K · Mₙ = '+format(sample.Mn_g_mol/1000,4)+' kg/mol';
    const inRange=(value,[a,b])=>value>=a&&value<=b;
    const isExt=sample.id===external?.id;
    const inSample=!!sample.T_range_K&&inRange(T,sample.T_range_K)&&inRange(point.phi,sample.phi_range);
    const osmTemps=(report?.osmotic?.points || []).map(row=>row.T_K).filter(Number.isFinite);
    const experimentalT=osmTemps.length?[Math.min(...osmTemps),Math.max(...osmTemps)]:null;
    const vdata=sys.calibration_range;
    $('domain').textContent='Bulutlanma kalibrasyonu '+format(vdata.T_K[0],5)+'–'+format(vdata.T_K[1],5)+' K ve φ = '+format(vdata.phi[0])+'–'+format(vdata.phi[1])+'. '+(isExt?'Fraksiyon II bağımsız basınç numunesidir.':inSample?'Bu seçim numunenin kayıtlı sıcaklık ve bileşim aralıklarında.':'Bu seçimin sıcaklığı veya bileşimi, numuneye ait bulutlanma kayıt aralığının dışında.')+($('model').value==='joint'&&experimentalT?' Ek ozmotik veri '+format(experimentalT[0],5)+'–'+format(experimentalT[1],5)+' K aralığında.':'');
    const a2Text=s.A2>0?'Seyreltik limitte A₂ pozitiftir ve net iki-zincir etkileşimi ideal katkıya göre basıncı artırır.':s.A2<0?'Seyreltik limitte A₂ negatiftir ve net iki-zincir etkileşimi ideal katkıya göre basıncı düşürür.':'İkinci virial terimi bu modelde sıfırdır.';
    $('regime').textContent=point.phase.count===2?'Bu ortalama derişimde dengedeki iki fazın bileşimleri φ₁ = '+scientific(s.eq.bin[0])+' ve φ₂ = '+scientific(s.eq.bin[1])+'. Basınç aynı çözücü kimyasal potansiyelinden hesaplanır. '+a2Text:a2Text;
    const params=new URLSearchParams({system:'ps',sample:sample.id,T:numberInput(T),phi:numberInput(point.phi)});
    const fitCode=$('model').value==='fh0'?'fh0':$('model').value==='cloud'?'fh1':'combined';
    params.set('fit',fitCode);
    const canTransfer=!isExt && $('model').value!=='joint' && point.phi>0 && T>=278 && T<=301;
    $('phaseLink').href=canTransfer?'phase.html?'+params:'phase.html?system=ps';
    $('phaseLink').textContent=canTransfer?'Faz diyagramında incele':'PS referans numunelerinin faz diyagramı';
    plot(s,point,scanLimit);
  }
  for (const id of inputIds) $(id).addEventListener(id==='sample'||id==='model'?'change':'input',()=>update(id));
  $('axis').addEventListener('change',()=>{if(latest)plot(latest.s,latest.point,scanLimit);});
  $('scan').addEventListener('input',()=>{$('concentration').value=numberInput(+$('scan').value); update('scan');});
  $('csv').addEventListener('click',()=>{
    if(!latest)return;
    const {s,sample,model} = latest;
    const headers=['sample_id','model','T_K','Mn_g_mol','Mw_g_mol','vbar_ml_g','Vs_ml_mol','integral_chi_phi0','integral_chi_slope','A2_mL_mol_g2','c_g_L','phi_polymer','phase_count','Pi_equilibrium_kPa','Pi_ideal_kPa','Pi_over_c_kPa_L_g','selected'];
    const lines=chartRows.map(row=>[sample.id,model,s.T,s.Mn,sample.Mw_g_mol,s.v,s.Vs,s.chi0,s.C,s.A2,row.c,row.phi,row.phase.count,row.pressure,row.ideal,row.reduced,Number(row.c===latest.point.c)].join(','));
    const url=URL.createObjectURL(new Blob([[headers.join(','),...lines].join('\n')+'\n'],{type:'text/csv;charset=utf-8;'}));
    const a=document.createElement('a');a.href=url;a.download='PS_siklohekzan_ozmotik_izoterm_'+numberInput(s.T)+'K.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
  $('volume').textContent='Derişim dönüşümünde φ = c v̄ / 1000 kullanılır (c g/L, v̄ = '+format(v,5)+' mL/g). Siklohekzan kafes referans hacmi '+format(Vs,6)+' mL/mol alınır ve sıcaklıkla sabit tutulur. Mₙ = Mw/Đ ve N = Mₙv̄/Vs.';
  let comparison=[];
  function updateExternal() {
    if(!external || !models.has('joint')) { $('external').hidden=true;return; }
    const selected=$('externalT').value;
    comparison=external.observations.filter(r=>selected==='all'||r.T_K===Number(selected)).map(r=>{
      const predictions=['cloud','joint'].map(key=>{
        const s=Osmosis.state(r.T_K,external.Mn_g_mol,v,Vs,models.get(key));
        return Osmosis.evaluate(r.c_g_L,s);
      });
      return {...r,cloud:predictions[0],joint:predictions[1]};
    });
    const meanError=key=>comparison.reduce((sum,r)=>sum+Math.abs(r[key].pressure-r.Pi_kPa),0)/comparison.length;
    $('score').textContent=comparison.length+' dış sınama ölçümünde ortalama mutlak hata: yalnız bulutlanma '+format(meanError('cloud'),3)+' kPa, bulutlanma + A₂ '+format(meanError('joint'),3)+' kPa.';
    const body=$('externalRows');body.replaceChildren();
    for(const r of comparison){
      const tr=body.insertRow();
      for(const cell of [format(r.T_K,5),format(r.c_g_L),format(r.Pi_kPa,4),format(r.cloud.pressure,4),format(r.joint.pressure,4),r.cloud.phase.count+' / '+r.joint.phase.count])tr.insertCell().textContent=cell;
    }
  }
  $('externalT').onchange=updateExternal;
  $('showSeries').onclick=()=>{
    const T=$('externalT').value==='all'?307.15:Number($('externalT').value);
    $('sample').value=external.id; $('temperature').value=T;
    $('concentration').value=external.observations.find(r=>r.T_K===T).c_g_L;
    update('sample'); $('plot').scrollIntoView({block:'center',behavior:'smooth'});
  };
  $('externalCsv').onclick=()=>{
    const heading='source_doi,page,table,fraction,Mn_g_mol,T_K,c_g_L,Q_1e6,Pi_measured_kPa,Pi_cloud_kPa,cloud_phase_count,Pi_joint_kPa,joint_phase_count';
    const lines=comparison.map(r=>[external.sourceDOI,external.sourcePage,external.sourceTable,'II',external.Mn_g_mol,r.T_K,r.c_g_L,r.Q1e6,r.Pi_kPa,r.cloud.pressure,r.cloud.phase.count,r.joint.pressure,r.joint.phase.count].join(','));
    const url=URL.createObjectURL(new Blob([[heading,...lines].join('\n')+'\n'],{type:'text/csv;charset=utf-8;'}));
    const a=document.createElement('a');a.href=url;a.download='PS_siklohekzan_Fraksiyon_II_basinc.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  updateExternal();update('initial');
})();
