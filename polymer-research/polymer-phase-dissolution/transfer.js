(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const R = window.TRANSFER_RESULTS;
  if (!R) { $('missing').hidden = false; return; }
  const finite = x => typeof x === 'number' && Number.isFinite(x);
  const fmt = (x, max = 3, min = 0) => {
    if (!finite(x)) return '—';
    if (x !== 0 && Math.abs(x) < 0.5 * 10 ** (-max)) {
      return x.toLocaleString('tr-TR', {notation:'scientific', maximumSignificantDigits:3});
    }
    return x.toLocaleString('tr-TR', {minimumFractionDigits:min, maximumFractionDigits:max});
  };
  const el = (tag, text, className) => {
    const node = document.createElement(tag);
    if (text !== undefined) node.textContent = text;
    if (className) node.className = className;
    return node;
  };
  const svg = (tag, attrs = {}) => {
    const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [key,val] of Object.entries(attrs)) node.setAttribute(key, val);
    return node;
  };
  const link = (text, href) => {
    const a = el('a',text);
    a.href = href;
    return a;
  };
  const markerDefs = {
    obs: {label:'Ölçüm',color:'#183a30'},
    cloud: {label:'Bulutlanma modeli',color:'#a35332'},
    combined: {label:'Bulutlanma + ozmotik veri',color:'#16685a'}
  };
  const modelIds = ['cloud','combined'];
  const tabText = {
    coex:['Ölçülen ve hesaplanan faz bileşimleri','Her sıcaklıktaki polimerce fakir ve zengin fazı kendi bileşim ekseninde karşılaştırın. Her ölçüm çifti uyum dışında tutulmuştur.'],
    osmotic:['İkinci ozmotik virial katsayısı','Seyreltik çözeltide polimer zincirleri arasındaki etkin etkileşim A₂ ile izlenir. Aynı sıcaklıkta farklı polimer fraksiyonlarının değerlerini seçebilirsiniz.'],
    cloud:['Bulutlanma sıcaklığı','İki modelin bulutlanma tahminlerini aynı numune ve bileşimde ölçülen sıcaklıkla karşılaştırın. Kalibrasyon ve dış sınama kayıtları ayrı gösterilir.']
  };
  const co = R.coexistence || {}, os = R.osmotic || {};
  let view = new URLSearchParams(location.search).get('view') || 'coex';
  if (!tabText[view]) view = 'coex';
  let rows = [], selected = null;
  const T = k => finite(k) ? ($('tempUnit').value === 'K' ? k : k-273.15) : null;
  const tUnit = () => $('tempUnit').value === 'K' ? 'K' : '°C';
  const temp = k => finite(k) ? fmt(T(k),3) + ' ' + tUnit() : '—';
  const cloudSources = [
    {title:'UChicago 3PDB · Bulutlanma ölçümleri',url:'https://pppdb.uchicago.edu/cloud_points'},
    {title:'Siporska ve ark. (2003) · Kalibrasyon',url:'https://doi.org/10.1039/B301190A'},
    {title:'Saeki ve ark. (1973) · Dış sınama',url:'https://doi.org/10.1021/ma60032a020'}
  ];
  const sourceSet = () => view === 'coex' ? co : view === 'osmotic' ? os : {sources:cloudSources,points:R.cloudPoints || []};
  const sourceTitle = s => s.title || s.label || s.name || s.doi || s.url || '';
  function addOption(select,value,text) {
    const opt = el('option',text); opt.value = value; select.append(opt);
  }
  function labelMass(raw) {
    const m = raw.Mw_g_mol || raw.Mn_g_mol;
    return finite(m) ? (raw.Mw_g_mol ? 'Mw' : 'Mₙ') + ' ' + fmt(m/1000,2) + ' kg/mol' : '';
  }
  function normalizedRows() {
    const set = sourceSet();
    const records = view==='osmotic'
      ? [...(set.points || []).map(row=>({...row,split:'calibration'})),
         ...(set.holdoutPoints || []).map(row=>({...row,split:'holdout'}))]
      : set.points || [];
    return records.map((raw,id) => {
      const sample = String(raw.series || raw.sample || raw.fraction || raw.Mw_g_mol || raw.Mn_g_mol || 'Numune');
      const sampleName = raw.sample_label || (view === 'osmotic' ? sample : labelMass(raw) || sample);
      const stamp = raw.split==='holdout'?' · dış sınama':raw.split==='calibration'?' · kalibrasyon':'';
      return {raw,id,sample,sampleName:sampleName+stamp};
    });
  }
  function activeModels() { return $('model').value === 'all' ? modelIds : [$('model').value]; }
  function pointValue(item,observable,model) {
    const row = item.raw;
    if (view === 'coex') return model === 'obs' ? row[observable] : row.predictions?.[model]?.[observable];
    if (view === 'osmotic') {
      const v = model === 'obs' ? row.A2_mol_ml_g2 : row.predictions?.[model];
      return finite(v) ? 10000*v : null;
    }
    return model === 'obs' ? row.T_K : row.predictions?.[model];
  }
  function asRecord(item) {
    const v = item.raw;
    return (item.sampleName || item.sample) + ' · ' + temp(v.T_K) + (view === 'cloud' && finite(v.phi_p) ? ' · φ '+fmt(v.phi_p,4) : '');
  }
  function drawSymbol(s,kind,x,y,r) {
    const def=markerDefs[kind];
    let shape;
    if (kind === 'obs') shape=svg('circle',{cx:x,cy:y,r,fill:'#fff',stroke:def.color,'stroke-width':1.6});
    else if (kind==='cloud') shape=svg('path',{d:'M'+x+' '+(y-r)+'l'+r+' '+r+'l-'+r+' '+r+'l-'+r+' -'+r+'z',fill:def.color});
    else shape=svg('rect',{x:x-r,y:y-r,width:2*r,height:2*r,rx:.55,fill:def.color});
    s.append(shape);
  }
  function ticks(min,max,num) {
    const rough=(max-min)/(num-1), decade=10**Math.floor(Math.log10(Math.max(rough,1e-12)));
    const step=[1,2,2.5,5,10].map(x=>x*decade).find(x=>x>=rough*.95) || 10*decade;
    const out=[];
    for(let i=Math.ceil((min-1e-11)/step);i*step<=max+1e-10;i++) out.push(Math.abs(i*step)<1e-11?0:i*step);
    return out;
  }
  function span(values,positive=false) {
    const lo=Math.min(...values),hi=Math.max(...values);
    let pad=(hi-lo)*.08 || Math.max(Math.abs(lo)*.015,.001);
    return [positive?Math.max(0,lo-pad):lo-pad,hi+pad];
  }
  function makeGraph(observable,title,xLabel,yLabel) {
    const kinds = view === 'cloud' ? activeModels() : ['obs',...activeModels()];
    const points = [];
    for(const item of rows) {
      const x = view === 'cloud' ? item.raw.phi_p : T(item.raw.T_K);
      for(const kind of kinds) {
        let y=pointValue(item,observable,kind);
        if (view==='cloud' && finite(y)) y-=item.raw.T_K;
        if (finite(x) && finite(y)) points.push({x,y,item,kind});
      }
    }
    if (!points.length) return;
    const w = view === 'coex' ? 332 : 610, h = view === 'coex' ? 250 : 270;
    const margins={left:view==='coex'?64:70,right:17,top:14,bottom:49};
    const [xmin,xmax]=span(points.map(p=>p.x),view==='cloud');
    const ys=points.map(p=>p.y);
    if(view === 'cloud' || view === 'osmotic') ys.push(0);
    const [ymin,ymax]=span(ys,view==='coex');
    const X=x=>margins.left+(x-xmin)/(xmax-xmin)*(w-margins.left-margins.right);
    const Y=y=>h-margins.bottom-(y-ymin)/(ymax-ymin)*(h-margins.top-margins.bottom);
    const graph=el('div',undefined,'graph');
    if (title) graph.append(el('h4',title));
    const chart=svg('svg',{viewBox:'0 0 '+w+' '+h,role:'img','aria-label':(title || $('plotTitle').textContent)+'. '+xLabel+' ve '+yLabel});
    const path=(x1,y1,x2,y2,stroke,dash)=>{const line=svg('line',{x1,y1,x2,y2,stroke}); if(dash)line.setAttribute('stroke-dasharray',dash);chart.append(line)};
    const label=(x,y,value,anchor='middle')=>{const n=svg('text',{x,y,'text-anchor':anchor,'font-family':'system-ui,sans-serif','font-size':12,fill:'#3b5547'});n.textContent=value;chart.append(n);return n};
    for(const y of ticks(ymin,ymax,5)) {path(X(xmin),Y(y),X(xmax),Y(y),'#dce4dc'); label(X(xmin)-8,Y(y)+4,fmt(y,5),'end');}
    for(const x of ticks(xmin,xmax,view==='coex'?4:6)) {path(X(x),Y(ymin),X(x),Y(ymin)+4,'#8aa094');label(X(x),Y(ymin)+20,fmt(x,4));}
    path(X(xmin),Y(ymax),X(xmin),Y(ymin),'#718b7c');path(X(xmin),Y(ymin),X(xmax),Y(ymin),'#718b7c');
    if(ymin<0&&ymax>0)path(X(xmin),Y(0),X(xmax),Y(0),'#92a59a','4 3');
    label((X(xmin)+X(xmax))/2,h-4,xLabel);
    const cy=(Y(ymin)+Y(ymax))/2;const yy=label(14,cy,yLabel);yy.setAttribute('transform','rotate(-90 14 '+cy+')');
    // A selected observation is circled; identical model symbols may overlap exactly.
    for(const p of points) {
      const x=X(p.x),y=Y(p.y);
      if(p.item.id===selected) chart.append(svg('circle',{cx:x,cy:y,r:6.4,fill:'none',stroke:'#6c7d6b','stroke-width':1,'stroke-dasharray':'2 2'}));
      const g=svg('g',{style:'cursor:pointer'});
      const tooltip=svg('title');tooltip.textContent=markerDefs[p.kind].label+' · '+asRecord(p.item)+' · '+fmt(p.y,5);g.append(tooltip);
      drawSymbol(g,p.kind,x,y,p.kind==='obs'?3.55:2.7);
      g.append(svg('circle',{cx:x,cy:y,r:5.5,fill:'transparent'}));
      g.addEventListener('click',()=>{selected=p.item.id;updateReport()});
      chart.append(g);
    }
    graph.append(chart);$('plots').append(graph);
  }
  function valuesList(items){
    const dl=el('dl');
    for(const [k,v] of items){dl.append(el('dt',k),el('dd',v));}
    return dl;
  }
  function reading() {
    const box=$('readout'); box.replaceChildren();
    const item=rows.find(r=>r.id===selected);
    if(!item) {box.append(el('p','Bu görünümde seçilebilir bir ölçüm bulunmuyor.')); return;}
    const row=item.raw;
    box.append(el('h4','Deney'),valuesList([['Sıcaklık',temp(row.T_K)], ...(labelMass(row)?[['Molar kütle',labelMass(row)]]:[])]));
    if(view==='coex') {
      if(finite(row.phi_total))box.append(valuesList([['Başlangıç bileşimi φ',fmt(row.phi_total,5)],['Ölçülen kritik sıcaklık',temp(row.Tc_measured_K)]]));
      box.append(valuesList([['Fakir faz · φ',fmt(row.phi_lean,6)],['Zengin faz · φ',fmt(row.phi_rich,6)]]));
      for(const id of activeModels()){
        const prediction=row.predictions?.[id];box.append(el('h4',markerDefs[id].label));
        if(finite(prediction?.phi_lean)&&finite(prediction?.phi_rich)) box.append(valuesList([['Fakir faz · φ',fmt(prediction.phi_lean,6)],['Zengin faz · φ',fmt(prediction.phi_rich,6)]]));
        else box.append(el('p',prediction?.phases===1?'Bu sıcaklıkta tek faz öngörülüyor.':'Bu kayıt için eşfaz çifti hesaplanamadı.'));
      }
    } else if(view==='osmotic') {
      if (finite(row.reported_halfwidth_mol_ml_g2)) box.append(el('p','Kaynağın bildirdiği ölçüm aralığı ±'+fmt(row.reported_halfwidth_mol_ml_g2*1e4,4)+' (A₂ × 10⁴).'));
      box.append(el('h4','A₂ × 10⁴ · mol·cm³/g²'), valuesList(['obs',...activeModels()].map(k=>[markerDefs[k].label,fmt(pointValue(item,'A2',k),4)])));
    } else {
      box.append(valuesList([['Polimer hacim kesri φ',fmt(row.phi_p,6)], ...activeModels().map(k=>[markerDefs[k].label,temp(pointValue(item,'T',k))])]));
    }
    const source=row.source || row.source_doi || row.source_id;
    const refs=sourceSet().sources || [];
    const cited=refs.find(s=>s.id===source || s.url==='https://doi.org/'+source || s.doi===source);
    if(cited?.url) box.append(link(sourceTitle(cited),cited.url));
    else if(typeof source==='string'&&/^10\.\d{4,9}\//i.test(source)) box.append(link('Ölçümün kaynağı', 'https://doi.org/'+source));
    else if(refs.length===1&&refs[0].url) box.append(link(sourceTitle(refs[0]),refs[0].url));
  }
  function table(head,body) {
    const div=el('div',undefined,'tableWrap'), t=el('table'), tr=t.createTHead().insertRow();
    for(const title of head){const th=el('th',title);th.scope='col';tr.append(th)}
    const tb=t.createTBody();
    for(const values of body) {const r=tb.insertRow();for(const value of values) r.insertCell().textContent=value;}
    div.append(t);return div;
  }
  const mean=vs=>vs.length?vs.reduce((s,n)=>s+n,0)/vs.length:null;
  const rmse=vs=>vs.length?Math.sqrt(mean(vs.map(v=>v*v))):null;
  function scoring() {
    const box=$('scores');box.replaceChildren();$('finding').hidden=true;
    if(!rows.length){$('scoreNote').textContent='Uygun gözlem sayısı sıfır. Hata ölçütü hesaplanamadı.';return;}
    if(view==='coex'){
      const both=rows.filter(({raw})=>modelIds.every(k=>finite(raw.predictions?.[k]?.phi_lean)&&finite(raw.predictions?.[k]?.phi_rich))&&finite(raw.phi_lean)&&finite(raw.phi_rich));
      const vals=activeModels().map(k=>{
        const pred=rows.filter(({raw})=>finite(raw.predictions?.[k]?.phi_lean)&&finite(raw.predictions?.[k]?.phi_rich)).length;
        const lean=both.map(({raw})=>Math.abs(raw.predictions[k].phi_lean-raw.phi_lean));
        const rich=both.map(({raw})=>Math.abs(raw.predictions[k].phi_rich-raw.phi_rich));
        const gap=rows.map(({raw})=>{
          const predicted=raw.predictions?.[k];
          const modelGap=finite(predicted?.phi_lean)&&finite(predicted?.phi_rich)?predicted.phi_rich-predicted.phi_lean:predicted?.phases===1?0:null;
          return finite(modelGap)&&finite(raw.phi_lean)&&finite(raw.phi_rich)? Math.abs(modelGap-(raw.phi_rich-raw.phi_lean)):null;
        }).filter(finite);
        return [markerDefs[k].label,pred+' / '+rows.length,fmt(mean(gap),5),fmt(mean(lean),5),fmt(mean(rich),5)];
      });
      const gapFor=k=>mean(rows.map(({raw:r})=>{const p=r.predictions[k];if(!p || p.phases!==1 && (!finite(p.phi_rich)||!finite(p.phi_lean)))return null;return Math.abs((p.phases===1?0:p.phi_rich-p.phi_lean)-(r.phi_rich-r.phi_lean));}).filter(finite));
      const jointCoverage=rows.filter(({raw:r})=>r.predictions.combined?.phases===2).length;
      const baseCoverage=rows.filter(({raw:r})=>r.predictions.cloud?.phases===2).length;
      const sharedError=k=>rmse(both.flatMap(({raw:r})=>[r.predictions[k].phi_rich-r.phi_rich,r.predictions[k].phi_lean-r.phi_lean]));
      $('finding').hidden=false;
      $('finding').textContent=rows.length+' eşfaz ölçümünün '+baseCoverage+' tanesinde bulutlanma modeli, '+jointCoverage+' tanesinde birleşik model iki faz öngördü. Eşfaz aralığının ortalama mutlak hatası sırasıyla '+fmt(gapFor('cloud'),4)+' ve '+fmt(gapFor('combined'),4)+'.'+(both.length?' İkisinin de iki faz verdiği '+both.length+' sıcaklıkta, iki uç birlikte değerlendirildiğinde RMSE sırasıyla '+fmt(sharedError('cloud'),4)+' ve '+fmt(sharedError('combined'),4)+'.':'');
      box.append(table(['Model','İki faz öngördüğü çift','Eşfaz aralığı · MAE (Δφ)','Fakir faz · MAE (φ)','Zengin faz · MAE (φ)'],vals));
      $('scoreNote').textContent='Eşfaz aralığı Δφ = φ(zengin) − φ(fakir). Tek faz öngörüsünde Δφ = 0 alınır. İki uç için ayrı MAE yalnız her iki modelin de eşfaz çifti verdiği '+both.length+' ortak kayıtta hesaplanır. MAE, ortalama mutlak hatadır.';
    } else {
      const vals=activeModels().map(k=>{
        const errors=rows.map(p=>{const v=pointValue(p,'',k),o=pointValue(p,'','obs');return finite(v)&&finite(o)?v-o:null}).filter(finite);
        return [markerDefs[k].label,String(errors.length),fmt(rmse(errors),3),fmt(mean(errors.map(Math.abs)),3),fmt(mean(errors),3)];
      });
      const unit=view==='osmotic'?'A₂ × 10⁴':'K';
      box.append(table(['Model','Kayıt','RMSE ('+unit+')','MAE ('+unit+')','Sapma ('+unit+')'],vals));
      $('scoreNote').textContent=view==='osmotic'?'Sapma, tahmin eksi ölçümdür. A₂ değerleri tabloda 10⁴ ile ölçeklenmiştir. Birleşik modelin uyumunda kullanılan kayıtların hatası kalibrasyon uyumunu gösterir.':'Sapma, tahmin eksi ölçümdür. Kalibrasyon ve bağımsız yayın ayrımını yukarıdaki numune seçiminden inceleyebilirsiniz.';
    }
  }
  function method() {
    const box=$('methodDetail');box.replaceChildren();
    const descriptions={
      coex:'Aynı sıcaklıktaki fakir ve zengin faz bileşimleri, modelin ortak teğet koşuluyla öngördüğü iki uçla karşılaştırılır. Faz bileşimleri boyutsuz polimer hacim kesridir.',
      osmotic:'Kütle derişimi c ile Π/(RT) = c/Mₙ + A₂c² + … . A₂ ölçümü, serbest enerjinin seyreltik uçtaki davranışını sınırlar.',
      cloud:'Bulutlanma ölçümü belirli bir başlangıç bileşimi ve sıcaklıktaki geçişi bildirir. İki eşfazın kimyasal bileşimi, ayrı bir deneyle ölçülebilir.'
    };
    box.append(el('p',descriptions[view]));
    const specs=R.models||[];
    if(specs.length){
      box.append(el('h4','Etkileşim fonksiyonu χ(φ,T) = A + B/T + Cφ'),
        table(['Model','A','B (K)','C'],specs.filter(m=>modelIds.includes(m.id)).map(m=>[markerDefs[m.id].label,fmt(m.A,5),fmt(m.B,4),fmt(m.C,5)])));
    }
    if(R.method && R.constants) {
      const c=R.constants;
      box.append(el('p','Birleşik uyumda iki veri grubu ayrı ortalanır. Sıcaklık kalıntıları '+fmt(c.scale_T_K,3)+' K, A₂ kalıntıları '+fmt(c.scale_A2_mol_ml_g2*1e4,2)+' × 10⁻⁴ mol·cm³/g² ile ölçeklenir. Bunlar uyum için seçilmiş ölçeklerdir.'),
        el('p','A₂ = (v̄²/Vₛ)[½ − A − B/T + C]. Referans kısmi özgül hacim v̄ = '+fmt(c.vbar_ml_g,4)+' cm³/g ve çözücünün molar hacmi Vₛ = '+fmt(c.V_s_ml_mol,3)+' cm³/mol.'));
    }
    if(view==='osmotic' && os.thetaHoldout) {
      const ref=os.thetaHoldout;
      box.append(el('h4','Teta sıcaklığının bağımsız kontrolü'),
        el('p','Işık saçılması ölçümünde A₂’nin sıfır olduğu sıcaklık '+temp(ref.T_K)+' (±'+fmt(ref.reported_halfwidth_K,2)+' K).'),
        table(['Model','Hesaplanan Θ','Ölçümden fark (K)'],specs.filter(m=>modelIds.includes(m.id)).map(m=>[markerDefs[m.id].label,temp(m.theta_K),fmt(m.thetaDelta_K,2)])),
        link('Bu modelden ozmotik basıncı hesapla','osmotic.html'));
    }
    if(view==='coex' && selected!==null){
      const row=rows.find(r=>r.id===selected)?.raw;
      if(row && finite(row.Tc_measured_K)){
        box.append(el('h4','Seçili numunenin kritik sıcaklığı'),table(['Kaynak','Kritik sıcaklık'],
          [['Ölçüm',temp(row.Tc_measured_K)],...specs.filter(m=>modelIds.includes(m.id)).map(m=>[markerDefs[m.id].label,temp(m.criticalBySeries?.[row.series]?.Tc_K)])]));
      }
    }
    if(view==='coex' && (R.sensitivity || []).some(r=>r.endpointMetrics?.n_total>0)) {
      box.append(el('h4','Ozmotik verinin ağırlığı değişirse'),
        table(['Ozmotik grubun ağırlığı','Bulutlanma dış RMSE (K)','Eşfaz aralığı MAE (Δφ)'],R.sensitivity.map(m=>[fmt(m.weight,2),fmt(m.cloudHoldout?.rmse_K,3),fmt(m.endpointMetrics?.mae_gap,5)])),
        el('p','Bu duyarlılık tablosu tüm dış sınama kayıtlarını kapsar. Ağırlık 0 yalnız bulutlanma uyumudur.'));
    }
    const set=sourceSet(), sources=set.sources || [];
    for(const note of [set.note,set.reason,set.limits])if(typeof note==='string'&&note)box.append(el('p',note));
    if(sources.length){box.append(el('h4','Ölçüm kaynakları'));const ul=el('ul');for(const s of sources){const li=el('li');if(typeof s==='string')li.textContent=s;else {const title=sourceTitle(s);if(s.url && /^https?:\/\//.test(s.url))li.append(link(title,s.url));else li.textContent=title;if(s.detail) li.append(document.createTextNode(' · '+s.detail));}ul.append(li);}box.append(ul);}
  }
  function recordChoices() {
    const o=$('observation');o.replaceChildren();
    for(const r of rows)addOption(o,String(r.id),asRecord(r));
    if(!rows.some(r=>r.id===selected)) selected=(view==='coex'?rows.at(-1):rows[0])?.id ?? null;
    o.value=selected===null?'':String(selected);o.disabled=!rows.length;
  }
  function updateReport() {
    const bank=normalizedRows(), choice=$('sample').value;
    rows=bank.filter(r=>choice==='*' || choice===r.sample || (choice==='@holdout'&&r.raw.split==='holdout') || (choice==='@calibration'&&r.raw.split==='calibration'));
    rows.sort((a,b)=>a.raw.T_K-b.raw.T_K || a.id-b.id);
    recordChoices();
    $('plots').replaceChildren();
    $('plotTitle').textContent={coex:'Polimer hacim kesri',osmotic:'Ozmotik virial katsayısı',cloud:'Bulutlanma sıcaklığındaki hata'}[view];
    $('figNote').textContent={coex:'Her simge bir sıcaklıktaki tek bir faz bileşimidir. İki grafiğin dikey ölçekleri farklıdır. Bir ölçümü seçmek için kaydı veya simgeyi kullanın.',osmotic:'Farklı numuneler aynı sıcaklıkta ölçülebilir. Aynı değeri alan simgeler üst üste gelir.',cloud:'Yatay sıfır çizgisi gözlenen bulutlanma sıcaklığıyla tam uyumu gösterir.'}[view];
    if (rows.length) {
      const xLabel='Sıcaklık ('+tUnit()+')';
      if(view==='coex')for(const [q,t] of [['phi_lean','Polimerce fakir faz'],['phi_rich','Polimerce zengin faz']])makeGraph(q,t,xLabel,'Hacim kesri φ');
      else if(view==='osmotic')makeGraph('A2','',xLabel,'A₂ × 10⁴ (mol·cm³/g²)');
      else makeGraph('T','','Polimer hacim kesri φ','Tahmin − ölçüm (K)');
    }
    $('plotEmpty').hidden=!!$('plots').childElementCount;
    $('plotEmpty').textContent='Bu seçime uygun karşılaştırılabilir ölçüm bulunmuyor.';
    $('legend').replaceChildren();
    for(const k of view==='cloud'?activeModels():['obs',...activeModels()]){const wrap=el('span');wrap.append(el('i',undefined,'swatch '+k),document.createTextNode(markerDefs[k].label));$('legend').append(wrap);}
    reading();scoring();method();
    $('download').disabled=!rows.length;
  }
  function choose(v) {
    view=v;selected=null;
    for(const b of document.querySelectorAll('[data-view]'))b.setAttribute('aria-pressed',String(b.dataset.view===v));
    [$('viewTitle').textContent,$('viewSub').textContent]=tabText[v];
    const sample=$('sample');sample.replaceChildren();addOption(sample,'*','Tümü');
    const bank=normalizedRows();
    if(view!=='coex'){for(const [id,text] of [['holdout','Bağımsız yayın'],['calibration','Kalibrasyon yayını']])if(bank.some(r=>r.raw.split===id))addOption(sample,'@'+id,text);}
    const groups=new Map(bank.map(r=>[r.sample,r.sampleName]));
    for(const [key,name] of groups)addOption(sample,key,name);
    sample.value=view!=='coex'&&bank.some(r=>r.raw.split==='holdout')?'@holdout':'*';
    sample.disabled=groups.size<2&&view==='coex';
    updateReport();
  }
  function save(name,text,type) {
    const href=URL.createObjectURL(new Blob([text],{type}));
    const a=document.createElement('a');a.href=href;a.download=name;a.click();
    setTimeout(()=>URL.revokeObjectURL(href),1000);
  }
  function toCSV() {
    const cols=['record','sample','T_K','Mw_g_mol','Mn_g_mol','split','observable','unit','polymer_volume_fraction','measurement','cloud_prediction','combined_prediction','source_reference'];
    const data=[cols];
    for(const item of rows){const r=item.raw;for(const observable of (view==='coex'?['phi_lean','phi_rich']:[view==='osmotic'?'A2':'Tcloud'])){
      const vals=['obs',...modelIds].map(k=>{const q=pointValue(item,observable,k);return !finite(q)?null:view==='osmotic'?q/10000:q});
      data.push([item.id,item.sample,r.T_K,r.Mw_g_mol,r.Mn_g_mol,r.split || (view==='coex'?'holdout':null),observable,view==='coex'?'volume fraction':view==='osmotic'?'mol cm3 g-2':'K',r.phi_p,...vals,r.source || r.source_id || (sourceSet().sources || []).map(s=>s.url || s.doi || s).join(' | ')]);
    }}
    const str=data.map(row=>row.map(v=>v==null?'':'"'+String(v).replaceAll('"','""')+'"').join(',')).join('\r\n')+'\r\n';
    save('PS_siklohekzan_'+view+'.csv',str,'text/csv;charset=utf-8');
  }
  $('download').addEventListener('click',toCSV);
  $('downloadAll').addEventListener('click',()=>save('PS_siklohekzan_model_kaydi.json',JSON.stringify(R,null,2),'application/json'));
  document.querySelectorAll('[data-view]').forEach(b=>b.addEventListener('click',()=>choose(b.dataset.view)));
  $('sample').addEventListener('change',()=>{selected=null;updateReport()});
  $('tempUnit').addEventListener('change',updateReport);$('model').addEventListener('change',updateReport);
  $('observation').addEventListener('change',()=>{selected=Number($('observation').value);updateReport()});
  $('study').hidden=false;choose(view);
})();
