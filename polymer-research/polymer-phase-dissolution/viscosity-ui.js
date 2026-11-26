(() => {
  'use strict';
  const $=id=>document.getElementById(id);
  const phys=window.Viscosity, catalog=window.VIS_INTRINSIC, external=window.VIS_EXTERNAL, finite=window.VIS_FINITE;
  if(!phys||!catalog?.rows?.length){$('error').textContent='Hesap motoru veya viskozite ölçüm tablosu yüklenemedi. Çevrimdışı paketi bütünüyle aç.';$('error').hidden=false;return;}
  const fmt=(x,n=4)=>Number.isFinite(x)?x.toLocaleString('tr-TR',{maximumSignificantDigits:n}):'—';
  const ns='http://www.w3.org/2000/svg';
  const opt=(s,k,label)=>s.add(new Option(label,k));
  const add=(p,t,attrs,txt)=>{const el=document.createElementNS(ns,t);for(const [k,v] of Object.entries(attrs))el.setAttribute(k,String(v));if(txt!=null)el.textContent=txt;p.append(el);return el;};
  const sourceLink=(parent,url,title)=>{if(!url)return;const a=document.createElement('a');a.href=url;a.textContent=title;parent.append(a);};
  const parse=(el)=>el.value.trim()===''?NaN:el.valueAsNumber;
  const csv=(name,lines)=>{
    const value=x=>'"'+String(x??'').replaceAll('"','""')+'"';
    const url=URL.createObjectURL(new Blob([lines.map(r=>r.map(value).join(',')).join('\n')+'\n'],{type:'text/csv;charset=utf-8'}));
    const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1200);
  };
  const q=new URLSearchParams(location.search);
  const queryNumber=k=>q.has(k)&&q.get(k).trim()!==''&&Number.isFinite(Number(q.get(k)))?Number(q.get(k)):null;
  const spec={
    cyclohexane:{id:'ps-cyclohexane-theta',label:'Siklohekzan',T:307.65},
    benzene:{id:'ps-benzene-25',label:'Benzen',T:298.15}
  };
  let hold=null, selection=null, chosenFinite=null;
  const context=()=>spec[$('solvent').value];
  const sampleRows=()=>catalog.rows.filter(r=>r.system===context().id&&r.study==='einaga1979');
  const refRows=()=>sampleRows().filter(r=>r.training);
  const externalRows=()=> (external?.rows||[]).filter(r=>r.system===context().id);
  function populate(prefer){
    const list=sampleRows();$('grade').replaceChildren();
    for(const r of list)opt($('grade'),r.id,r.id+' · '+fmt(r.Mw_g_mol/1000,5)+' kg/mol');
    const picked=list.find(r=>r.id===prefer)||list.find(r=>r.id==='BK2500-4')||list.find(r=>r.matched_Rg&&Number.isFinite(r.kH))||list[0];
    $('grade').value=picked.id;
    if(!$('mass').value && !q.has('source'))$('mass').value=String(picked.Mw_g_mol/1000);
  }
  const T=queryNumber('T'),Mw=queryNumber('Mw'),c=queryNumber('c');
  const solvent=q.get('solvent');
  if(spec[solvent])$('solvent').value=solvent;
  populate(q.get('sample'));
  if(c!==null)$('concentration').value=String(c);
  let transferredTwoPhase=q.get('phaseCount')==='2';
  const handed=['source','sample','system','Mw','Mn','T','solvent'].some(k=>q.has(k));
  if(handed){
    const matched=sampleRows().find(r=>r.id===q.get('sample'));
    const validSystem=q.get('system')==='ps' || (!q.has('system')&&!!matched&&!!spec[solvent]);
    const validT=T!==null && Math.abs(T-context().T)<1e-5;
    const identity=validSystem&&spec[solvent]&&validT;
    const massAgrees=Mw===null||[matched?.Mw_g_mol,matched?.Rg_Mw_g_mol].some(m=>Number.isFinite(m)&&Math.abs(m-Mw)<=1e-7*m);
    const gradeMatch=identity&&matched&&massAgrees;
    if(gradeMatch)$('mode').value='measured';
    else if(identity && Mw>0)$('mode').value='mh';
    else $('mode').value='custom';
    if(Mw>0)$('mass').value=String(Mw/1000);
    const parts=[q.get('sample')||q.get('name'),solvent?spec[solvent]?.label||solvent:null,T!==null?fmt(T,5)+' K':null,c!==null?fmt(c)+' g/L':null].filter(Boolean);
    $('handoff').textContent='Aktarılan koşul · '+(parts.join(' · ')||'kullanıcı numunesi')+'. '+(gradeMatch?'Viskozite ölçümü aynı fraksiyon ve çözücü koşuluyla eşleşti.':$('mode').value==='mh'?'Bu fraksiyon için eşlenmiş viskozite ölçümü bulunmuyor. Tahmin seçtiğin koşuldaki mol kütlesi uyumuna dayanır.':'Bu koşul için eşlenmiş bir viskozite ölçümü bulunmuyor. Kendi numunene ait [η] ölçümünü girebilirsin.');
    $('handoff').hidden=false;
  }
  function digits(values,n){
    const out=values.map(v=>fmt(v,n));
    return new Set(out).size===out.length?out:values.map(v=>fmt(v,5));
  }
  function tickValues(min,max,log) {
    if(log){
      const all=[];
      const first=Math.floor(Math.log10(min)),last=Math.ceil(Math.log10(max));
      for(let i=first;i<=last;i++)for(const v of [1,2,5]){const x=v*10**i;if(x>=min&&x<=max) all.push(x);}
      if(all.length>=3)return all.filter((v,i)=>i%Math.ceil(all.length/6)===0);
    }
    return Array.from({length:5},(_,i)=>min+i*(max-min)/4);
  }
  function plot(svg,xmin,xmax,ymin,ymax,xName,yName,log=false) {
    svg.replaceChildren();const l=75,r=588,t=28,b=243;
    const v=x=>log?Math.log(x):x;
    const X=x=>l+(r-l)*(v(x)-v(xmin))/(v(xmax)-v(xmin)),Y=y=>b-(b-t)*(v(y)-v(ymin))/(v(ymax)-v(ymin));
    const text=(x,y,s,props={})=>add(svg,'text',{x,y,fill:'#445e51','font-size':'14','font-family':'system-ui,sans-serif',...props},s);
    const ys=tickValues(ymin,ymax,log), ylabels=digits(ys,3);
    ys.forEach((y,i)=>{add(svg,'line',{x1:l,x2:r,y1:Y(y),y2:Y(y),stroke:'#dae2d9'});text(l-8,Y(y)+4,ylabels[i],{'text-anchor':'end'});});
    add(svg,'path',{d:'M'+l+' '+t+' V'+b+' H'+r,fill:'none',stroke:'#91aa96'});
    const xs=tickValues(xmin,xmax,log),labels=digits(xs,3);
    xs.forEach((x,i)=>text(X(x),b+23,labels[i],{'text-anchor':'middle'}));
    text((l+r)/2,291,xName,{'text-anchor':'middle'});text(l,t-11,yName,{'font-weight':600});
    function path(data,color,dash='none') {
      const d=data.filter(p=>p.every(Number.isFinite)).map(([x,y],i)=>(i?'L':'M')+X(x).toFixed(2)+' '+Y(y).toFixed(2)).join(' ');
      add(svg,'path',{d,fill:'none',stroke:color,'stroke-width':2.5,'stroke-dasharray':dash});
    }
    function dot(x,y,kind,label,selected=false) {
      if(!Number.isFinite(x)||!Number.isFinite(y))return;
      if(selected)add(svg,'circle',{cx:X(x),cy:Y(y),r:9,fill:'none',stroke:'#966433','stroke-width':1.6});
      const item=kind==='additional'?add(svg,'circle',{cx:X(x),cy:Y(y),r:3.9,fill:'#a87832',stroke:'white','stroke-width':1}):kind==='external'?add(svg,'rect',{x:X(x)-4,y:Y(y)-4,width:8,height:8,fill:'white',stroke:'#9b5930','stroke-width':2,transform:'rotate(45 '+X(x)+' '+Y(y)+')'}):add(svg,'circle',{cx:X(x),cy:Y(y),r:4.3,fill:'#273f32',stroke:'white','stroke-width':1});
      add(item,'title',{},label);
    }
    return {path,dot};
  }
  function renderSample() {
    const mode=$('mode').value, custom=mode==='custom', measured=mode==='measured', system=context();
    $('solventField').hidden=custom;$('gradeField').hidden=!measured;$('massField').hidden=mode!=='mh';
    $('customEtaField').hidden=!custom;$('customKField').hidden=!custom;
    const grade=measured?sampleRows().find(r=>r.id===$('grade').value):null,mw=measured?grade.Mw_g_mol:custom?(Mw>0?Mw:null):parse($('mass'))*1000;
    const concentration=parse($('concentration'));
    const training=refRows();const fitted=phys.fitMH(training);
    const eta=measured?grade.eta_mL_g:custom?parse($('customEta')):(mw>0&&Number.isFinite(mw)?phys.predictMH(mw,fitted):null);
    const kH=grade?grade.kH:custom?parse($('customK')):null;
    const validation=mode==='mh' && !(mw>0&&Number.isFinite(mw))?'M𝑤 için kg/mol cinsinden pozitif bir değer gir.':!(eta>0&&Number.isFinite(eta))?'Numunene ait [η] değerini mL/g cinsinden gir.':!(concentration>=0&&Number.isFinite(concentration))?'Derişim için sıfır veya pozitif, sonlu bir değer gir.':custom&&$('customK').value.trim()&&!Number.isFinite(kH)?'Ölçülen kH değerini gir veya alanı boş bırak.':'';
    $('invalid').hidden=!validation;$('invalid').textContent=validation;$('sampleResults').hidden=!!validation;$('sampleCsv').disabled=!!validation;selection=null;
    const back=new URLSearchParams({source:'viscosity',system:'ps',T:String(system.T)});
    if(concentration>=0&&Number.isFinite(concentration))back.set('c',String(concentration));
    if(measured&&grade.matched_Rg){back.set('sample',grade.id);back.set('Mw',String(mw));}
    if(mode==='mh'&&mw>0)back.set('Mw',String(mw));
    if(custom){back.set('system','custom');if(T===null)back.delete('T');else back.set('T',String(T));}
    $('transfer').hidden=!(measured&&grade.matched_Rg);
    if(!$('transfer').hidden)$('transfer').href='viscosity-transfer.html?'+new URLSearchParams({sample:grade.id,solvent:$('solvent').value,...(Number.isFinite(concentration)&&concentration>=0?{c:String(concentration)}:{})});
    $('chainReturn').href=!custom&&$('solvent').value==='benzene'&&grade?.matched_Rg?'chain-solvents.html?'+new URLSearchParams({sample:grade.id,c:String(concentration)}):'chain.html?'+back;
    const aiParams=new URLSearchParams({source:'viscosity',system:custom?'custom':'ps'});
    if(!custom){aiParams.set('T',String(system.T));aiParams.set('solvent',$('solvent').value);if(grade)aiParams.set('sample',grade.id);}
    if(concentration>=0&&Number.isFinite(concentration))aiParams.set('c',String(concentration));
    $('ai').href='viscosity-ai.html'+(aiParams.size?'?'+aiParams:'');
    if(validation)return;
    const {u,c_eta_g_L}=phys.scale(eta,concentration);
    if(!Number.isFinite(u)||!Number.isFinite(c_eta_g_L)||c_eta_g_L<=0){$('invalid').textContent='Bu değer sayısal hesap aralığının dışında.';$('invalid').hidden=false;$('sampleResults').hidden=true;$('sampleCsv').disabled=true;return;}
    const kKnown=Number.isFinite(kH);
    let prediction=null;
    if(u<1&&kKnown){try{prediction=phys.dilute(eta,kH,concentration);}catch(e){if(!(e instanceof RangeError))throw e;}}
    let rel=concentration===0?1:prediction?.eta_rel??null;
    if(!Number.isFinite(rel)||rel<=0)rel=null;
    const outside=mode==='mh'&&(mw<fitted.minMw||mw>fitted.maxMw);
    $('eta').textContent=fmt(eta)+' mL/g';
    $('etaFrom').textContent=measured?'Seyreltik limit · '+grade.id+(Number.isFinite(grade.eta_pm_mL_g)?' · bildirilen ±'+fmt(grade.eta_pm_mL_g)+' mL/g':''):custom?'Kullanıcının ölçümü':outside?'Mark–Houwink · mol kütlesi uyum aralığının dışında':'Mark–Houwink · mol kütlesinden tahmin';
    $('u').textContent=fmt(u);$('ceta').textContent='cη = 1/[η] = '+fmt(c_eta_g_L)+' g/L';
    $('relative').classList.toggle('textual',rel===null);$('relative').textContent=rel!==null?fmt(rel):u>=1?'Seyreltik sınır aşıldı':!kKnown?'kH ölçümü gerekli':'Model kapsamı dışında';
    $('relativeLabel').textContent=concentration===0?'Saf çözücü sınırı':rel===null?'Sonlu derişimde η/η₀':'Huggins açılımı · kH = '+fmt(kH);
    $('interpretation').textContent=concentration===0?'Seçili noktada η = η₀. [η], polimer derişimi sıfıra yaklaşırken bağıl viskozitenin başlangıç eğimini verir.':u>=1?'Bu ekranda [η]c ≥ 1 için ikinci mertebe açılımından sayısal bağıl viskozite verilmez. Bu gösterim sınırı, akışta ani bir geçişi veya dolaşıklığı tanımlamaz.':!kKnown?'Bu numune ve koşul için bir Huggins katsayısı tanımlanmamış. İntrinsik viskozite ve viskozimetrik derişim ölçeği hesaplanabilir.':rel===null?'Girilen katsayı ve derişim fiziksel bir bağıl viskozite vermiyor.':'Bu değer ikinci mertebeye kadar seyreltik çözelti öngörüsüdür. [η]c küçüldükçe yüksek mertebeli katkılar azalır. Derişime özgü deneysel hata sınırı verilmemiştir.';
    if(transferredTwoPhase&&c!==null&&Math.abs(concentration-c)<1e-9){
      rel=null;$('relative').textContent='—';$('relativeLabel').textContent='İki fazlı toplam bileşim';$('interpretation').textContent='Kaynak denge modeli iki faz öngörüyor. Hacimsel ortalama derişim için tek bir homojen sıvı viskozitesi hesaplanamaz. Ayrı fazların derişimleriyle incele.';
    }
    const pure=!custom?pureViscosity($('solvent').value,system.T):null;
    showDynamic($('dynamic'),pure,rel);
    selection={mode,grade,mw,concentration,eta,u,c_eta_g_L,kH,rel,T:custom?T:system.T,pure};
  }
  function pureViscosity(solvent,T_K){
    return window.VIS_SOLVENT?.records?.find(r=>r.solvent===solvent&&Math.abs(r.T_K-T_K)<1e-6&&r.p_Pa===101325);
  }
  function showDynamic(el,pure,rel){
    el.hidden=!(pure&&rel>0&&Number.isFinite(pure.eta_mPa_s));el.replaceChildren();
    if(el.hidden)return;
    el.append('Dinamik viskozite η ≈ '+fmt(rel*pure.eta_mPa_s)+' mPa·s. Saf çözücü için η₀ = '+fmt(pure.eta_mPa_s)+' mPa·s, '+fmt(pure.T_K,5)+' K ve 101,325 kPa. ');
    sourceLink(el,pure.source.url,pure.source.label);
  }
  function tableRows(id,lines){
    const parent=$(id);parent.replaceChildren();
    for(const cells of lines){
      const row=document.createElement('tr');
      for(const s of cells){const td=document.createElement('td');td.textContent=s;row.append(td);}
      parent.append(row);
    }
  }
  function renderValidation(){
    const train=refRows(), outside=externalRows(), sys=context();
    const fit=phys.fitMH(train);
    const report=phys.studyBenchmark(train,outside);
    const separate=$('split').value==='external'&&!!report.external;
    hold=separate?report.external:report.loo;
    const shown=separate?outside:[];
    const name=id=>shown.find(r=>r.id===id)?.sample||id;
    const metric=hold.metrics, folds=hold.folds, low=fit.minMw/1000, high=fit.maxMw/1000;
    const extrap=folds.filter(r=>r.domain==='extrapolation').length;
    $('exponent').textContent=fmt(fit.a);$('rms').textContent='%'+fmt(metric.rmse_pct,3);$('mae').textContent=fmt(metric.mae_mL_g)+' mL/g';
    $('mhSolvent').textContent='Polistiren · '+sys.label+' · '+fmt(sys.T-273.15)+' °C';
    $('splitNote').textContent=(separate?'Dış yayından '+folds.length+' ölçüm. ':'Her adımda bir fraksiyon çıkarılarak model yeniden uyduruldu. ')+extrap+'/'+folds.length+' sınama noktası ilgili uyumun kütle aralığının dışında. Ana uyum '+fmt(low)+'–'+fmt(high)+' kg/mol.';
    $('extKey').hidden=!separate;
    const selectedExtra=selection?.grade&&!selection.grade.training?selection.grade:null;
    $('extraKey').hidden=!selectedExtra;
    const points=[...train,...shown,...(selectedExtra?[selectedExtra]:[])];
    let xmin=Math.min(...points.map(r=>r.Mw_g_mol/1000))*.86, xmax=Math.max(...points.map(r=>r.Mw_g_mol/1000))*1.16;
    const predicted=x=>phys.predictMH(1000*x,fit);
    const values=[...points.map(r=>r.eta_mL_g),predicted(xmin),predicted(xmax)];
    const chart=plot($('mhPlot'),xmin,xmax,Math.min(...values)*.83,Math.max(...values)*1.18,'M𝑤 (kg/mol) · log','[η] (mL/g) · log',true);
    chart.path([[xmin,predicted(xmin)],[low,predicted(low)]],'#246e51','5 4');
    chart.path([[low,predicted(low)],[high,predicted(high)]],'#246e51');
    chart.path([[high,predicted(high)],[xmax,predicted(xmax)]],'#246e51','5 4');
    for(const row of points){
      const x=row.Mw_g_mol/1000, y=row.eta_mL_g;
      chart.dot(x,y,row===selectedExtra?'additional':shown.includes(row)?'external':'train',name(row.id)+' · '+fmt(x)+' kg/mol · '+fmt(y)+' mL/g',selection?.grade?.id===row.id);
    }
    tableRows('fitRows',folds.map(r=>[name(r.id),fmt(r.Mw_g_mol/1000,5),fmt(r.measured_mL_g),fmt(r.predicted_mL_g),r.domain==='interpolation'?'İçinde':'Dışında']));
    const refs=$('intrinsicSources');refs.replaceChildren();refs.append('Uyum için ');
    sourceLink(refs,catalog.source.url,catalog.source.name+' · Tablo I');
    if(separate){
      for(const src of external.sources.filter(s=>shown.some(r=>r.study===s.id))){refs.append(' · Sınama için ');sourceLink(refs,src.url,src.title);}
    }
    refs.append('. Grafik tüm uyum fraksiyonlarından elde edilen eğriyi gösterir. Kesikli uzantılar mol kütlesi aralığının dışındadır.');
  }
  function finiteY(rel,c,kind,intercept){
    if(kind==='rel')return rel;
    if(c===0)return intercept;
    return (kind==='inh'?Math.log(rel):rel-1)/(c/1000);
  }
  function renderFinite(){
    const series=finite?.series?.find(s=>s.id===$('series').value);
    if(!series){$('finiteError').textContent='Kapiler ölçüm tablosu yüklenemedi.';$('finiteError').hidden=false;$('finiteResults').hidden=true;return;}
    const all=phys.finiteBenchmark(series.rows),split=$('finiteSplit').value;
    const report=split==='highest'?all:all[split];
    if(!report){$('finiteError').textContent='Bu sınama için yeterli sayıda farklı derişim bulunmuyor.';$('finiteError').hidden=false;$('finiteResults').hidden=true;return;}
    $('finiteError').hidden=true;$('finiteResults').hidden=false;chosenFinite={series,report};
    const {fit,folds,metrics}=report;
    const tests=folds.filter(r=>r.role==='test'),last=tests.at(-1);
    $('finiteCondition').textContent=series.label;
    $('finiteTitle').textContent='Kapiler ölçüm ve iki seyreltik çözelti modeli';
    $('heldTitle').textContent='Sınama · '+fmt(last.c_g_L)+' g/L';
    $('finiteObs').textContent=fmt(last.eta_rel);$('finitePred').textContent=fmt(last.huggins_eta_rel);$('finiteK').textContent=fmt(last.kraemer_eta_rel);
    const relation=split==='lowInterpolation'?'İki düşük derişimle uyum, aradaki derişimle sınama. İki parametreli uyumun artık serbestlik derecesi sıfır.':'Yüksek derişimde dışa uzatım. ';
    $('finiteNote').textContent=relation+' Sınanan '+tests.length+' noktada bağıl RMSE · Huggins %'+fmt(metrics.huggins.rmse_pct,3)+' · Kraemer %'+fmt(metrics.kraemer.rmse_pct,3)+'. Gösterilen noktada [η]c = '+fmt(last.u)+'.';
    const kind=$('ordinate').value;
    const observed=folds.map(r=>finiteY(r.eta_rel,r.c_g_L,kind));
    const maxC=Math.max(...folds.map(r=>r.c_g_L)),h=[],k=[];
    for(let i=0;i<=80;i++){
      const c=maxC*i/80, C=c/1000;
      const relH=1+C*(fit.eta_mL_g+fit.hugginsSlope_mL2_g2*C);
      const relK=Math.exp(C*(fit.kraemerIntercept_mL_g+fit.kraemerSlope_mL2_g2*C));
      h.push([c,finiteY(relH,c,kind,fit.eta_mL_g)]);
      k.push([c,finiteY(relK,c,kind,fit.kraemerIntercept_mL_g)]);
    }
    const ys=[...observed,...h.map(p=>p[1]),...k.map(p=>p[1])],yMin=Math.min(...ys),yMax=Math.max(...ys),pad=Math.max(1e-4,(yMax-yMin)*.1);
    const label={rel:'Bağıl viskozite ηᵣ',red:'ηₛₚ/c (mL/g)',inh:'ln(ηᵣ)/c (mL/g)'}[kind];
    const chart=plot($('finitePlot'),0,maxC*1.055,yMin-pad,yMax+pad,'Polimer derişimi (g/L)',label);
    chart.path(h,'#246e51');chart.path(k,'#89628b','6 4');
    for(const [i,row] of folds.entries())chart.dot(row.c_g_L,observed[i],row.role==='test'?'external':'train',fmt(row.c_g_L)+' g/L · '+label+' = '+fmt(observed[i]));
    $('finiteMethod').textContent='Huggins kesişimi [η] = '+fmt(fit.eta_mL_g)+' mL/g, kH = '+fmt(fit.kH)+'. Kraemer kesişimi [η] = '+fmt(fit.kraemerIntercept_mL_g)+' mL/g, işaretli kK = '+fmt(fit.kK)+'. '+series.measurementNote;
    const lines=series.rows.map(row=>{
      const calc=folds.find(r=>r.id===row.id);
      return [fmt(row.c_g_L),fmt(row.eta_rel),fmt(calc?.huggins_eta_rel),fmt(calc?.kraemer_eta_rel),fmt(row.solutionViscosity_mPa_s),!calc?'Kullanılmadı':calc.role==='test'?'Sınama':'Uyum'];
    });
    tableRows('finiteRows',lines);
    $('finiteAi').href='viscosity-ai.html?'+new URLSearchParams({source:'viscosity',series:series.id,c:String(last.c_g_L),T:String(series.T_K),system:'ps'});
    const fs=$('finiteSource');fs.replaceChildren();
    sourceLink(fs,series.source.url,series.source.label);
    if(series.molarMassNote)fs.append(' · '+series.molarMassNote);
    const fluid=series.solvent.toLocaleLowerCase('tr-TR');
    const pure=pureViscosity(fluid==='toluen'?'toluene':fluid.includes('siklohekzan')?'cyclohexane':fluid,series.T_K);
    showDynamic($('finiteDynamic'),pure,last.huggins_eta_rel);
  }
  const selectView=showFinite=>{
    $('sampleView').hidden=showFinite;$('finiteView').hidden=!showFinite;
    $('sampleTab').setAttribute('aria-pressed',String(!showFinite));
    $('finiteTab').setAttribute('aria-pressed',String(showFinite));
    if(showFinite)renderFinite();
  };
  const refresh=()=>{renderSample();renderValidation();};
  for(const id of ['grade','mode','solvent','mass'])$(id).addEventListener('change',()=>{transferredTwoPhase=false;});
  $('solvent').addEventListener('change',()=>{populate($('grade').value);refresh();});
  for(const id of ['grade','mode'])$(id).addEventListener('change',refresh);
  for(const id of ['mass','concentration','customEta','customK'])$(id).addEventListener('input',renderSample);
  $('split').addEventListener('change',renderValidation);
  for(const id of ['series','ordinate','finiteSplit'])$(id).addEventListener('change',renderFinite);
  $('sampleTab').addEventListener('click',()=>selectView(false));$('finiteTab').addEventListener('click',()=>selectView(true));
  $('sampleCsv').addEventListener('click',()=>{
    if(!selection)return;
    const s=selection;
    csv('viskozite-numune.csv',[['sample','mode','solvent','T_K','Mw_g_mol','c_g_L','intrinsic_mL_g','reduced_concentration','c_eta_g_L','kH','eta_rel_Huggins','eta0_mPa_s_at_101325Pa','solution_mPa_s_model','source'],[
      s.grade?.id||'',s.mode,s.mode==='custom'?q.get('solvent'):$('solvent').value,s.T,s.mw,s.concentration,s.eta,s.u,s.c_eta_g_L,Number.isFinite(s.kH)?s.kH:null,s.rel,s.pure?.eta_mPa_s,s.pure&&s.rel!==null?s.pure.eta_mPa_s*s.rel:null,s.mode==='custom'?'user':catalog.source.doi]]);
  });
  $('intrinsicCsv').addEventListener('click',()=>{
    if(!hold)return;
    csv('intrinsik-viskozite-sinama.csv',[['id','study','T_K','solvent','Mw_g_mol','observed_mL_g','predicted_mL_g','domain','training_ids','fit_min_Mw_g_mol','fit_max_Mw_g_mol'],...hold.folds.map(r=>[r.id,r.study,context().T,$('solvent').value,r.Mw_g_mol,r.measured_mL_g,r.predicted_mL_g,r.domain,r.train_ids.join('|'),r.training_min_Mw,r.training_max_Mw])]);
  });
  $('finiteCsv').addEventListener('click',()=>{
    if(!chosenFinite)return;
    const {series,report}=chosenFinite;
    csv('kapiler-viskozite-sinama.csv',[['series','T_K','id','c_g_L','solution_time_s','solvent_time_s','reported_dynamic_mPa_s','eta_rel_reported','eta_rel_Huggins','eta_rel_Kraemer','reduced_concentration','role','source_url','reported_concern'],...series.rows.map(row=>{
      const r=report.folds.find(f=>f.id===row.id);
      return [series.id,series.T_K,row.id,row.c_g_L,row.solutionTime_s,row.solventTime_s,row.solutionViscosity_mPa_s,row.eta_rel,r?.huggins_eta_rel,r?.kraemer_eta_rel,r?.u,r?.role||'unused',series.source.url,row.reportedConcern===true];
    })]);
  });
  for(const series of finite?.series||[])opt($('series'),series.id,series.label);
  refresh();selectView(q.get('view')==='finite');$('app').hidden=false;
})();
