(() => {
  'use strict';
  const el=id=>document.getElementById(id);
  const fmt=(n,d=4)=>Number.isFinite(n)?n.toLocaleString('tr-TR',{maximumSignificantDigits:d}):'—';
  const gp=window.ViscosityGP, series=window.VIS_FINITE?.series, intrinsic=window.VIS_INTRINSIC?.rows;
  if(!gp||!series?.length||!intrinsic){el('error').textContent='Viskozite modelleri veya kaynak ölçümler yüklenemedi.';el('error').hidden=false;return;}
  const report=gp.buildBenchmark(series,intrinsic);
  const protocols={wholeGrade:'Seçili mol kütlesinin tüm ölçümleri',highest:'En yüksek derişim',upperTwo:'En yüksek iki derişim',lowInterpolation:'Alt aralıktaki bir ara derişim'};
  const q=new URLSearchParams(location.search);
  const start=series.find(s=>s.id===(q.get('series')||q.get('sample')));
  const defaultSample=start||series.find(s=>s.id.endsWith('633000'))||series[0];
  for(const sample of series)el('series').add(new Option(sample.label,sample.id));
  el('series').value=defaultSample.id;
  let sample,record,fit;
  const percent=(x,y)=>100*(x/y-1);
  function probe(){
    const c=el('probe').valueAsNumber, limit=Math.max(...sample.rows.map(r=>r.c_g_L));
    if(!(Number.isFinite(c)&&c>=0&&c<=limit)){
      ['probeBase','probeGP','probeObs'].forEach(id=>el(id).textContent='—');
      el('support').textContent='Bu deney karşılaştırması için 0–'+fmt(limit)+' g/L arasında bir derişim gir.';return;
    }
    const p=gp.predict(fit,c), observed=sample.rows.find(r=>Math.abs(r.c_g_L-c)<1e-8);
    el('probeBase').textContent=fmt(p.baseline);el('probeGP').textContent=fmt(p.gp);el('probeObs').textContent=fmt(observed?.eta_rel);
    el('support').textContent=c===0?'Saf çözücü sınırı her iki modelde tam olarak 1.':
      'u = '+fmt(p.u)+'. '+(p.domain==='interpolation'?'Eğitimdeki boyutsuz derişim aralığında.':'Eğitimdeki boyutsuz derişim aralığının dışında.')+
      (p.u>=1?' u ≥ 1; seyreltik yaklaşım burada bir ekstrapolasyon.':'')+
      (fit.etaPriorDomain==='extrapolation'?' [η] için mol kütlesi de uyum aralığının dışında.':'')+
      (!p.physicallyAdmissible?' Model bu noktada kabul edilen pozitif ve artan dal koşulunu karşılamıyor.':'');
  }
  function chart(){
    const svg=el('chart');svg.replaceChildren();
    const add=(name,attrs,txt)=>{const e=document.createElementNS('http://www.w3.org/2000/svg',name);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,String(v));if(txt!=null)e.textContent=txt;svg.append(e);return e;};
    const same=record.strategy!=='wholeGrade';
    const xhi=Math.max(...record.folds.map(r=>r.c_g_L),...(same?fit.trainingObservations.map(r=>r.c_g_L):[]));
    const points=Array.from({length:101},(_,i)=>gp.predict(fit,i*xhi/100));
    const yhi=1+1.13*(Math.max(1.1,...points.flatMap(p=>[p.gp,p.baseline]),...record.folds.map(r=>r.observed),...(same?fit.trainingObservations.map(r=>r.observed):[]))-1);
    const X=x=>61+462*x/xhi,Y=y=>237-196*(y-1)/(yhi-1);
    const txt=(x,y,t,attrs={})=>add('text',{x,y,fill:'#485c54','font-size':14,...attrs},t);
    for(let i=0;i<=4;i++){const y=1+(yhi-1)*i/4;add('line',{x1:61,x2:523,y1:Y(y),y2:Y(y),stroke:'#dfe6de'});txt(54,Y(y)+5,fmt(y,3),{'text-anchor':'end'});txt(X(xhi*i/4),260,fmt(xhi*i/4,3),{'text-anchor':'middle'});}
    add('path',{d:'M61 33V237H523',fill:'none',stroke:'#859a8f'});
    txt(62,20,'Bağıl viskozite ηᵣ');txt(302,285,'Polimer derişimi (g/L)',{'text-anchor':'middle'});
    for(const [k,color,dash] of [['baseline','#19694f','none'],['gp','#915b36','7 5']]){
      const d=points.map((p,i)=>(i?'L':'M')+X(p.c_g_L).toFixed(2)+' '+Y(p[k]).toFixed(2)).join(' ');
      add('path',{d,fill:'none',stroke:color,'stroke-width':2.6,'stroke-dasharray':dash});
    }
    if(same) for(const p of fit.trainingObservations)add('circle',{cx:X(p.c_g_L),cy:Y(p.observed),r:3.8,fill:'#244b39'});
    for(const p of record.folds){const x=X(p.c_g_L),y=Y(p.observed);add('path',{d:'M'+x+' '+(y-5)+'l5 5 -5 5 -5 -5z',stroke:'#92532c',fill:'white','stroke-width':2});}
    const key=svg.nextElementSibling.querySelector('.pt');key.parentElement.hidden=!same;
  }
  function display(keepProbe=false){
    sample=series.find(s=>s.id===el('series').value);
    const strategy=el('protocol').value;
    record=[...report.within,...report.transfer].find(r=>r.testSeriesId===sample.id&&r.strategy===strategy);
    const train=series.find(s=>s.id===record.trainSeriesId),whole=strategy==='wholeGrade';
    fit=whole?gp.fitTransfer(train,sample,intrinsic):gp.fitWithin(sample,strategy);
    el('baseScore').textContent=fmt(record.metrics.baseline.rmse);el('gpScore').textContent=fmt(record.metrics.gp.rmse);
    el('count').textContent=record.folds.length+' nokta';el('countNote').textContent=whole?'Bir numunenin tüm serisi':'Bu deneyden ayrılan derişim';
    el('condition').textContent=sample.label;
    el('scope').textContent=whole?'Grafikte hedef numunenin ölçümleri gösterilir. Diğer numunede öğrenilen katsayılar bu numuneye aktarılır. İki numune aynı tezde aynı çözücü ve sıcaklıkta ölçülmüştür.':sample.measurementNote;
    el('splitInfo').textContent='Eğitim serisi · '+train.label+' ('+record.trainIds.length+' nokta). [η] = '+fmt(fit.eta_mL_g)+' mL/g. '+fit.etaPriorSource+'.'+(strategy==='lowInterpolation'?' İki eğitim noktası Huggins doğrusunu tam belirler. Gauss sürecinin öğrenebileceği artık bu bölmede sıfırdır.':'');
    el('source').textContent=sample.source.label;el('source').href=sample.source.url;
    const body=el('held');body.replaceChildren();
    for(const f of record.folds){
      const row=body.insertRow();for(const v of [fmt(f.c_g_L),fmt(f.observed),fmt(f.baseline),fmt(f.gp),fmt(percent(f.baseline,f.observed),3)+' %',fmt(percent(f.gp,f.observed),3)+' %'])row.insertCell().textContent=v;
    }
    if(!keepProbe)el('probe').value=String(record.folds.at(-1).c_g_L);
    el('probe').max=String(Math.max(...sample.rows.map(r=>r.c_g_L)));probe();chart();
  }
  function options(){
    const old=el('protocol').value;
    const eligible=report.transfer.some(r=>r.testSeriesId===el('series').value);
    const sel=el('protocol');sel.replaceChildren();
    for(const [k,v] of Object.entries(protocols))if(k!=='wholeGrade'||eligible)sel.add(new Option(v,k));
    if(Array.from(sel.options).some(o=>o.value===old))sel.value=old;
    display();
  }
  el('download').onclick=()=>{
    const header='sample,split,training_sample,training_rows,held_row,c_g_L,u,observed_eta_rel,Huggins_eta_rel,GP_eta_rel,conditional_latent_GP_SD,concentration_domain,eta_prior_domain';
    const lines=record.folds.map(r=>[sample.id,record.strategy,fit.trainSeriesId,fit.trainIds.join('|'),r.id,r.c_g_L,r.u,r.observed,r.baseline,r.gp,r.gpConditionalSD,r.domain,fit.etaPriorDomain].join(','));
    const url=URL.createObjectURL(new Blob([[header,...lines].join('\n')+'\n'],{type:'text/csv;charset=utf-8'}));
    const a=document.createElement('a');a.href=url;a.download='viskozite_AI_'+sample.id+'_'+record.strategy+'.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  function show(){el('unmatched').hidden=true;el('experiment').hidden=false;options();}
  el('useReference').onclick=show;el('series').onchange=options;el('protocol').onchange=()=>display();el('probe').oninput=probe;
  const wantsMatch=['source','system','series','sample','solvent','T','Mw','name'].some(k=>q.has(k));
  const solventOK=!q.has('solvent')||((start?.solvent==='Toluen'?'toluene':'cyclohexane')===q.get('solvent'));
  const match=start&&(!q.has('system')||q.get('system')==='ps')&&solventOK&&(!q.has('T')||Math.abs(+q.get('T')-start.T_K)<1e-5);
  if(wantsMatch&&!match){
    el('unmatched').hidden=false;el('unmatchedText').textContent='Aktarılan numune veya deney koşulu için sonlu derişimde eşlenmiş AI eğitim serisi bulunmuyor. Referans deneylerden birini ayrı olarak açabilirsin.';
  }else{show();if(match&&q.has('c')&&q.get('c').trim()){el('probe').value=q.get('c');probe();}}
})();
