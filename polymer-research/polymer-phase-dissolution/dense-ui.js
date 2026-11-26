(function () {
  'use strict';

  const $ = id => document.getElementById(id);
  const results = window.DENSE_RESULTS;
  const fh = window.DenseFH;
  const ns = 'http://www.w3.org/2000/svg';
  const palette = {base:'#685581', extended:'#146950', guard:'#2e6f8d', observed:'#a74628', grid:'#dce4dc', paper:'#fffefa'};
  const pretty = (n, digits=2) => Number.isFinite(n)
    ? n.toLocaleString('tr-TR', {maximumFractionDigits:digits, minimumFractionDigits:digits}) : '—';
  const concentration = n => n.toLocaleString('tr-TR',{maximumFractionDigits:2});
  const celsius = t => (Math.round((t-273.15)*100)/100).toLocaleString('tr-TR');
  const getFold = () => results.folds.find(item => item.id === $('fold').value);
  const getTemp = () => Number($('temp').value);
  const series = () => getFold().rows.filter(r=>r.T_K===getTemp()).sort((a,b)=>a.c_g_L-b.c_g_L);
  const getRecord = () => series().find(r=>String(r.id)===$('record').value);
  const auxiliary = row => row.source_tier === 'secondary_concentration';
  const sampleRole = row => auxiliary(row) ? 'Sınama · ikincil derişim' : row.role === 'test' ? 'Sınama' : 'Uyum';
  const guardShown = () => $('guard').checked;
  const variants = () => guardShown() ? ['base','extended','guard'] : ['base','extended'];
  function shortVirial(x) {
    return x === 0 ? '0' : x.toLocaleString('tr-TR',{notation:'scientific',maximumFractionDigits:2});
  }
  let view = 'pressure';
  let chartStates = new Map();

  function element(name, attrs, body) {
    const node = document.createElementNS(ns,name);
    for (const [key,value] of Object.entries(attrs||{})) node.setAttribute(key,String(value));
    if (body !== undefined) node.textContent=body;
    return node;
  }
  function options(select, items) {
    select.replaceChildren(...items.map(([value,label])=>{
      const option=document.createElement('option');
      option.value=value;
      option.textContent=label;
      return option;
    }));
  }
  function error(message) {
    $('dataError').textContent=message;
    $('dataError').hidden=false;
    $('workspace').hidden=true;
  }
  if (!results?.folds?.length || !fh?.state || !fh?.evaluate) {
    error('Karşılaştırma verileri veya denge çözücüsü yüklenemedi. Çevrimdışı kullanımda paketin tamamını aynı dizine çıkarın.');
    return;
  }

  function onFold() {
    const fold=getFold();
    chartStates.clear();
    const temperatures=[...new Set(fold.rows.map(r=>r.T_K))].sort((a,b)=>a-b);
    const previousT=getTemp();
    options($('temp'),temperatures.map(t=>[t,celsius(t)+' °C']));
    const requested=Number(new URLSearchParams(location.search).get('T'));
    const activeTemps=[...new Set(fold.rows.filter(r=>r.role==='test').map(r=>r.T_K))];
    const target=activeTemps.includes(requested)?requested:activeTemps.includes(previousT)?previousT:activeTemps.at(-1);
    if (temperatures.includes(target)) $('temp').value=String(target);
    $('tier').value='primary';
    $('tierChoice').hidden=fold.score_ids.length===fold.test_ids.length;
    renderScore();
    const formula='A₃ = (⅓ − 2C + 2D) v̄³/Vₛ';
    $('virialText').textContent='Ek terimli model '+formula+' ile '+shortVirial(fold.A3)+' mol mL²/g³ verir. Bağımsız PS fraksiyonlarında θ yakınında (≈34,5 °C) ışık saçılması A₃ için pozitif değerler gösterir. İşaret kısıtı uygulanınca uyum A₃ = 0 sınırında kalıyor. Basınç hatası ve seyreltik çözeltideki bu bağımsız gözlem birlikte değerlendirilmelidir.';
    $('fitText').textContent='D, uyuma ayrılan mutlak basınçların kPa cinsinden karesel sapmasını en aza indirir. Arama aralığı −1 ≤ D ≤ 1/12. Üst sınır bu modelde en fazla bir eşfaz çifti sağlayan matematiksel bir seçimdir. A₃ ≥ 0 koşuluyla tekrar uydurulan katsayı D = '+pretty(fold.D_guard,5)+'.';
    const fit=results.base, s=results.sample;
    $('paramText').textContent='Seçili uyumda A = '+pretty(fit.A,4)+', B = '+pretty(fit.B,3)+' K, C = '+pretty(fit.C,4)+' ve D = '+pretty(fold.D,5)+'. Referans çözücü molar hacmi '+pretty(s.Vs_ml_mol,2)+' mL/mol, polimer özgül hacmi '+pretty(s.vbar_ml_g,3)+' mL/g.';
    $('sampleName').textContent='Fraksiyon II · Mₙ '+pretty(s.Mn_g_mol/1000,0)+' kg/mol';
    onTemperature();
  }
  function renderScore() {
    const f=getFold(), all=$('tier').value==='all';
    const scores=all ? f.metrics_all : f.metrics;
    const count=(all ? f.test_ids : f.score_ids).length;
    $('baseMae').textContent=pretty(scores.base.mae)+' kPa';
    $('newMae').textContent=pretty(scores.extended.mae)+' kPa';
    $('guardMae').textContent=pretty(scores.guard.mae)+' kPa';
    $('guardScore').hidden=!guardShown();
    $('guardCol').hidden=!guardShown();
    const labels=f.test_ids.length===f.score_ids.length
      ? 'D için '+f.train_ids.length+' uyum ölçümü. Hata uyum dışında bırakılan '+count+' ölçüm üzerinden hesaplandı.'
      : all ? 'İki ikincil derişim kaydı dahil üç sıcaklık-serisi ölçümünün duyarlılık sonucu.'
            : 'D için '+f.train_ids.length+' uyum ölçümü. İki kaynakta teyit edilen tek derişimdeki mutlak hata gösteriliyor. Diğer iki derişim duyarlılık hesabında kullanılabilir.';
    $('scoreText').textContent=labels;
  }
  function onTemperature() {
    const visible=series();
    options($('record'),visible.map(r=>[r.id,concentration(r.c_g_L)+' g/L']));
    const lastTest=visible.filter(r=>r.role==='test' && !auxiliary(r)).at(-1) || visible.filter(r=>r.role==='test').at(-1);
    $('record').value=String((lastTest||visible.at(-1)).id);
    $('tableNote').textContent=celsius(getTemp())+' °C izotermi · '+visible.length+' ölçüm. Bir derişim seçerek grafikteki karşılığını inceleyin.';
    renderRecord();
  }
  function field(name,value,css='') {
    const div=document.createElement('div'), dt=document.createElement('dt'), dd=document.createElement('dd');
    div.className=css; dt.textContent=name; dd.textContent=value; div.append(dt,dd); return div;
  }
  function renderRecord() {
    const r=getRecord();
    if (!r) return;
    $('sampleRole').textContent=auxiliary(r)
      ? 'Sınama ölçümü · derişim ikincil kaynaktan'
      : r.role==='test'?'Sınama ölçümü · D uyumunun dışında':'Uyum ölçümü · D değerinin belirlenmesinde';
    $('values').replaceChildren(
      field('Ölçülen Π',pretty(r.observed_kPa,3)+' kPa','observed'),
      field('Temel model',pretty(r.base_kPa,3)+' kPa'),
      field('Ek terimli model',pretty(r.extended_kPa,3)+' kPa'),
      ...(guardShown() ? [field('A₃ ≥ 0 modeli',pretty(r.guard_kPa,3)+' kPa')] : []),
      field('Polimer hacim kesri φ',pretty(r.c_g_L*results.sample.vbar_ml_g/1000,4))
    );
    $('phaseNote').textContent='Denge faz sayısı: temel model '+r.base_phases+', ek terimli model '+r.extended_phases+(guardShown()?', A₃ ≥ 0 modeli '+r.guard_phases:'')+'.';
    const body=$('observations');
    body.replaceChildren(...series().map(item=>{
      const row=document.createElement('tr');
      if (String(item.id)===String(r.id)) row.className='is-selected';
      const key=document.createElement('td'), action=document.createElement('button');
      action.textContent=concentration(item.c_g_L);
      action.setAttribute('aria-label',action.textContent+' g/L ölçümünü incele');
      action.setAttribute('aria-pressed',String(item.id===r.id));
      action.addEventListener('click',()=>selectRecord(item.id));
      key.append(action);row.append(key);
      for (const value of [sampleRole(item),pretty(item.observed_kPa,3),pretty(item.base_kPa,3),pretty(item.extended_kPa,3),...(guardShown()?[pretty(item.guard_kPa,3)]:[])]) {
        const td=document.createElement('td'); td.textContent=value; row.append(td);
      }
      return row;
    }));
    draw();
  }
  function selectRecord(id) {
    $('record').value=String(id);
    renderRecord();
  }
  function state(T,D) {
    const key=T+'|'+D;
    if (!chartStates.has(key)) {
      const s=results.sample;
      chartStates.set(key,fh.state(T,s.Mn_g_mol,s.vbar_ml_g,s.Vs_ml_mol,results.base,D));
    }
    return chartStates.get(key);
  }
  function legend(entries) {
    $('legend').replaceChildren(...entries.map(([style,label])=>{
      const key=document.createElement('span'), symbol=document.createElement('i');
      key.className='key'; symbol.className=style; symbol.setAttribute('aria-hidden','true');
      key.append(symbol,document.createTextNode(label));
      return key;
    }));
  }
  function tickStep(max,steps) {
    const power=10**Math.floor(Math.log10(Math.max(max,1e-8)/steps));
    const m=max/steps/power;
    return (m>5?10:m>2.5?5:m>2?2.5:m>1?2:1)*power;
  }
  function tickLabel(n) {
    return n.toLocaleString('tr-TR',{maximumFractionDigits:3});
  }
  function line(svg,attrs) {svg.append(element('line',attrs));}
  function draw() {
    const host=$('plot'), record=getRecord();
    if (!record || !host.clientWidth) return;
    const fold=getFold(), visible=series(), T=getTemp();
    const width=Math.max(236,Math.min(850,host.clientWidth)), compact=width<430;
    const height=compact?250:285;
    const dims={l:compact?39:53,r:width-14,t:28,b:height-49};
    const axisLabel=view==='pressure'?'Π (kPa)':'Tahmin − ölçüm (kPa)';
    const svg=element('svg',{viewBox:'0 0 '+width+' '+height,role:'img','aria-label':axisLabel+'; '+celsius(T)+' °C’de polimer derişiminin fonksiyonu.'});
    svg.append(element('title',{},view==='pressure'?'Ölçülen basınç ve iki denge modeli':'Her bir ölçümde işaretli model sapmaları'));
    const xDataMax=Math.max(...visible.map(r=>r.c_g_L));
    const xStep=tickStep(xDataMax,compact?4:6);
    const xMax=Math.ceil(xDataMax/xStep)*xStep;
    const x=c=>dims.l+c/xMax*(dims.r-dims.l);
    const pressure=view==='pressure';
    const curves={}, residuals=[];
    let yMin=0, yMax=0;
    if (pressure) {
      for (const key of variants()) {
        const st=state(T,{base:0,extended:fold.D,guard:fold.D_guard}[key]);
        curves[key]=Array.from({length:171},(_,i)=>{
          const c=xDataMax*i/170;
          const pi=fh.evaluate(c,st).pressure;
          if (!Number.isFinite(pi)) throw new RangeError('Denge basıncı sonlu olmalı.');
          yMax=Math.max(pi,yMax);
          return [c,pi];
        });
      }
      yMax=Math.max(yMax,...visible.map(r=>r.observed_kPa));
    } else {
      for (const r of visible) {
        for (const key of variants()) {
          const value=r[key+'_kPa']-r.observed_kPa;
          residuals.push({r,key,value});
          yMin=Math.min(yMin,value);yMax=Math.max(yMax,value);
        }
      }
    }
    const yStep=tickStep(Math.max(yMax-yMin,0.01),4);
    const yLo=Math.floor(yMin/yStep)*yStep, yHi=Math.max(yLo+yStep,Math.ceil(yMax/yStep)*yStep);
    const y=v=>dims.b-(v-yLo)/(yHi-yLo)*(dims.b-dims.t);
    const font={'font-size':compact?11:12};
    for (let t=yLo;t<=yHi+yStep/100;t+=yStep) {
      const tick=Math.abs(t)<yStep/1000?0:t;
      line(svg,{x1:dims.l,x2:dims.r,y1:y(tick),y2:y(tick),stroke:tick===0?'#93a798':palette.grid});
      svg.append(element('text',{x:dims.l-7,y:y(tick)+4,'text-anchor':'end',...font},tickLabel(tick)));
    }
    for (let t=0;t<=xMax+xStep/100;t+=xStep) {
      line(svg,{x1:x(t),x2:x(t),y1:dims.b,y2:dims.b+4,stroke:'#95aa99'});
      svg.append(element('text',{x:x(t),y:dims.b+18,'text-anchor':'middle',...font},tickLabel(t)));
    }
    svg.append(element('text',{x:4,y:14,'font-size':compact?12:13,'font-weight':650},axisLabel),
      element('text',{x:(dims.l+dims.r)/2,y:height-5,'text-anchor':'middle','font-size':compact?12:13},'Polimer derişimi, c (g/L)'));
    line(svg,{x1:x(record.c_g_L),x2:x(record.c_g_L),y1:dims.t,y2:dims.b,stroke:'#acb9ac','stroke-dasharray':'2 4'});
    if (pressure) {
      for (const key of variants()) {
        const path=curves[key].map(([c,p],i)=>(i?'L':'M')+x(c).toFixed(2)+' '+y(p).toFixed(2)).join(' ');
        svg.append(element('path',{d:path,fill:'none',stroke:palette[key],'stroke-width':2.5,'stroke-dasharray':key==='base'?'6 4':key==='guard'?'2 4':'none'}));
      }
      for (const r of visible) {
        const at=x(r.c_g_L), yp=y(r.observed_kPa);
        if (r.id===record.id) svg.append(element('circle',{cx:at,cy:yp,r:7.6,fill:'none',stroke:'#483b29','stroke-width':1.2}));
        const point=auxiliary(r)
          ? element('path',{d:'M'+at+' '+(yp-5.5)+'l5.5 5.5 -5.5 5.5 -5.5 -5.5Z',fill:palette.paper,stroke:palette.observed,'stroke-width':1.7})
          : element('circle',{cx:at,cy:yp,r:4.4,fill:r.role==='test'?palette.observed:palette.paper,stroke:palette.observed,'stroke-width':1.7});
        point.append(element('title',{},concentration(r.c_g_L)+' g/L · '+pretty(r.observed_kPa,3)+' kPa · '+sampleRole(r)));
        point.style.cursor='pointer';point.addEventListener('click',()=>selectRecord(r.id));svg.append(point);
      }
      legend([['base','Temel model'],['extended','Ek terimli model'],...(guardShown()?[['guard','A₃ ≥ 0 modeli']]:[]),
        ...['test','train'].filter(role=>visible.some(r=>r.role===role && !auxiliary(r))).map(role=>[role,role==='test'?'Sınama ölçümü':'Uyum ölçümü']),
        ...(visible.some(auxiliary)?[['aux','Sınama · ikincil derişim']]:[])]);
      $('plotNote').textContent='Eğriler seçili sıcaklıkta termodinamik dengeden hesaplanır. Halkayla çevrili nokta seçili membran osmometrisi ölçümüdür.';
    } else {
      for (const {r,key,value} of residuals) {
        const cx=x(r.c_g_L), cy=y(value);
        const shape=key==='base'?element('path',{d:'M'+cx+' '+(cy-5)+'l5 5 -5 5 -5 -5Z',fill:palette.paper,stroke:palette.base,'stroke-width':1.8})
          :key==='guard' ? element('rect',{x:cx-3.4,y:cy-3.4,width:6.8,height:6.8,fill:'none',stroke:palette.guard,'stroke-width':1.7})
          :element('circle',{cx,cy,r:3,fill:palette.extended,stroke:palette.extended});
        shape.append(element('title',{},concentration(r.c_g_L)+' g/L · '+pretty(value,3)+' kPa · '+sampleRole(r)));
        shape.style.cursor='pointer';shape.addEventListener('click',()=>selectRecord(r.id));svg.append(shape);
      }
      legend([['basePoint','Temel model'],['extPoint','Ek terimli model'],...(guardShown()?[['guardPoint','A₃ ≥ 0 modeli']]:[])]);
      $('plotNote').textContent='Pozitif sapma basıncın yüksek tahmin edildiğini gösterir. Uyuma ve sınamaya ayrılan ölçümler alttaki tabloda belirtilir.';
    }
    $('plotTitle').textContent=pressure?'Basınç izotermi':'Ölçüme göre sapma';
    $('plotSub').textContent=celsius(T)+' °C · '+visible.filter(r=>r.role==='test').length+' sınama, '+visible.filter(r=>r.role==='train').length+' uyum ölçümü';
    host.replaceChildren(svg);
  }

  function exportCsv() {
    const fold=getFold(), rows=series();
    const header=['fold','role','source_tier','primary_score','record','T_K','c_g_L','Mn_g_mol','D','D_guard','Pi_observed_kPa','Pi_base_kPa','Pi_extended_kPa','Pi_guard_kPa','resid_base_kPa','resid_extended_kPa','resid_guard_kPa','model_base_phases','model_extended_phases','model_guard_phases','source_DOI'];
    const data=rows.map(r=>[fold.id,r.role,r.source_tier,fold.score_ids.includes(r.id),r.id,r.T_K,r.c_g_L,results.sample.Mn_g_mol,fold.D,fold.D_guard,r.observed_kPa,r.base_kPa,r.extended_kPa,r.guard_kPa,r.base_kPa-r.observed_kPa,r.extended_kPa-r.observed_kPa,r.guard_kPa-r.observed_kPa,r.base_phases,r.extended_phases,r.guard_phases,'10.1021/ja01517a022']);
    const escape=v=>'"'+String(v).replaceAll('"','""')+'"';
    const csv=[header,...data].map(row=>row.map(escape).join(',')).join('\r\n')+'\r\n';
    const blob=new Blob([csv],{type:'text/csv;charset=utf-8;'});
    const url=URL.createObjectURL(blob), a=document.createElement('a');
    a.href=url;a.download='PS_Siklohekzan_'+fold.id+'_'+getTemp()+'K.csv';
    document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  $('pressureTab').addEventListener('click',()=>setView('pressure'));
  $('residualTab').addEventListener('click',()=>setView('residual'));
  function setView(next) {
    view=next;
    $('pressureTab').setAttribute('aria-pressed',String(next==='pressure'));
    $('residualTab').setAttribute('aria-pressed',String(next==='residual'));
    draw();
  }
  $('export').addEventListener('click',exportCsv);
  $('tier').addEventListener('change',renderScore);
  $('guard').addEventListener('change',()=>{renderScore();renderRecord();});
  $('fold').addEventListener('change',onFold);
  $('temp').addEventListener('change',onTemperature);
  $('record').addEventListener('change',renderRecord);
  options($('fold'),results.folds.map(f=>[f.id,f.label]));
  const requestedFold=new URLSearchParams(location.search).get('fold');
  const defaultFold=[requestedFold,'hold44','high'].find(id=>results.folds.some(f=>f.id===id));
  if (defaultFold) $('fold').value=defaultFold;
  $('workspace').hidden=false;
  onFold();
  if (window.ResizeObserver) new ResizeObserver(()=>draw()).observe($('plot'));
  else window.addEventListener('resize',draw);
})();
