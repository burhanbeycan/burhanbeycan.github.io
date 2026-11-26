(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const phys = window.ChainPhysics, catalog = window.CHAIN_DATA, other=window.CHAIN_EXTERNAL;
  if (!phys || !catalog?.rows?.length) {
    $('loadError').textContent='Zincir modeli veya ölçüm tablosu yüklenemedi. Çevrimdışı paketi bütünüyle açın.';
    $('loadError').hidden=false;
    return;
  }
  const rows=catalog.rows.filter(r=>r.role!=='external'), ext=other?.rows || [];
  const all=[...rows,...ext], free=phys.fit(rows), theta=phys.fit(rows,.5);
  const loo=phys.benchmark(rows), external=ext.length ? phys.external(rows,ext) : null;
  const F=(x,n=4)=>Number.isFinite(x)?x.toLocaleString('tr-TR',{maximumSignificantDigits:n}):'—';
  const option=(sel,v,txt)=>sel.add(new Option(txt,v));
  const query=new URLSearchParams(location.search);
  const maybe=k=>query.has(k)&&query.get(k).trim()!==''&&Number.isFinite(Number(query.get(k)))?Number(query.get(k)):null;
  const inRefT=T=>Math.abs(T-catalog.T_K)<1e-5;
  const scaleName = row => row.source_id===rows[0].source_id?'Işık saçılması':'SAXS';
  all.forEach(r=>option($('sample'),r.id,r.label+' · '+F(r.Mw_g_mol/1000,5)+' kg/mol · '+scaleName(r)));
  const matched = all.find(r=>r.id===query.get('sample'));
  $('sample').value=matched?.id||rows[Math.floor(rows.length/2)].id;
  const supplied=query.has('sample')||query.has('source')||query.has('system')||query.has('Mw')||query.has('Mn')||query.has('Rg');
  const knownMass=maybe('Mw'), radius=maybe('Rg'), T=maybe('T'), mn=maybe('Mn');
  const c=maybe('c'); if(c!==null)$('c').value=String(c);
  if (supplied) {
    if (matched && (!T||inRefT(T)) && (!query.has('system')||query.get('system')==='ps')) $('mode').value='measured';
    else if (knownMass>0 && !radius && query.get('system')==='ps' && T!==null && inRefT(T)) $('mode').value='model';
    else $('mode').value='custom';
    if(knownMass>0)$('mw').value=String(knownMass/1000);
    if(radius>0)$('rg').value=String(radius);
    const origin=query.get('sample')||query.get('name')||'seçili numune';
    const context=[T!==null?'T = '+F(T,5)+' K':'',mn>0?'Mₙ = '+F(mn/1000)+' kg/mol':'',knownMass>0?'M𝑤 = '+F(knownMass/1000)+' kg/mol':'',c!==null?'c = '+F(c)+' g/L':''].filter(Boolean);
    const missing=mn>0 && !(knownMass>0) && !matched?' Aktarılan Mₙ ile M𝑤 farklı mol kütlesi ortalamalarıdır. Örtüşme hesabı için M𝑤 gerekir.':'';
    const tempNote= T!==null&&!inRefT(T)?' Referans saçılma serisi 34,5 °C’dedir. Kendi numunenin çözücü ve sıcaklığına ait Rg ölçümünü kullan.':'';
    const isTwo=maybe('phaseCount')===2, lean=maybe('cLean'), rich=maybe('cRich');
    const phaseNote=isTwo?' Kaynak model iki faz öngörüyor. Aktarılan c toplam derişimdir. Her fazın örtüşme oranı kendi derişimine göre değerlendirilmelidir.'+(lean!==null&&rich!==null?' Kaynak modelde c(fakir) = '+F(lean)+' ve c(zengin) = '+F(rich)+' g/L.':''):'';
    $('transferNote').textContent='Aktarılan koşul · '+origin+(context.length?' · '+context.join(' · '):'')+'.'+missing+tempNote+phaseNote;
    $('transferNote').hidden=false;
  }
  if(!external)$('validation').options[1].disabled=true;
  const ns='http://www.w3.org/2000/svg';
  const add=(parent,tag,attrs,txt)=>{
    const el=document.createElementNS(ns,tag);
    for(const [k,v] of Object.entries(attrs))el.setAttribute(k,String(v));
    if(txt!=null)el.textContent=txt;
    parent.append(el);return el;
  };
  function exportCsv(name,rows) {
    const cell=x=>'"'+String(x??'').replaceAll('"','""')+'"';
    const text=rows.map(r=>r.map(cell).join(',')).join('\n')+'\n';
    const url=URL.createObjectURL(new Blob([text],{type:'text/csv;charset=utf-8'})), a=document.createElement('a');
    a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  function logAxisTicks(lo,hi,maxTicks=6) {
    const values=[];
    for(let k=Math.floor(Math.log10(lo));k<=Math.ceil(Math.log10(hi));k++) {
      for (const f of [1,2,5]) {const x=f*10**k;if(x>=lo&&x<=hi)values.push(x);}
    }
    const spacing=Math.ceil(values.length/maxTicks);
    return values.filter((x,i)=>i%spacing===0);
  }
  let state=null, result=null, seed=17, xyz=null, coordsKey='', timer=null;
  let yaw=.3,pitch=-.3,drag=null;
  function coilPaint() {
    if(!xyz||$('chainPanel').hidden||$('computed').hidden)return;
    const canvas=$('coil'),ctx=canvas.getContext('2d');if(!ctx)return;
    const width=canvas.clientWidth||630,height=canvas.clientHeight||260;
    const dpr=Math.min(window.devicePixelRatio||1,2);
    canvas.width=width*dpr;canvas.height=height*dpr;ctx.scale(dpr,dpr);
    const cs=Math.cos(yaw),ss=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch);
    const model=xyz.points.map(([x,y,z])=>{const w=x*ss+z*cs;return [x*cs-z*ss, y*cp-w*sp,y*sp+w*cp];});
    const rad=Math.max(xyz.targetRg_nm*2.45,...model.map(p=>Math.hypot(...p)));
    const scale=.44*Math.min(height,width)/rad;
    const pts=model.map(([x,y,z])=>[width*.5+x*scale,height*.48-y*scale,z]);
    const segments=pts.slice(1).map((p,i)=>[pts[i],p]).sort((a,b)=>(a[0][2]+a[1][2])-(b[0][2]+b[1][2]));
    ctx.lineCap='round';
    for(const [a,b] of segments) {
      const t=Math.min(1,Math.max(0,.5+(a[2]+b[2])/(4*rad)));
      ctx.strokeStyle='rgb('+Math.round(62+t*170)+' '+Math.round(128+t*110)+' '+Math.round(98+t*63)+')';
      ctx.lineWidth=2.1+t*2;ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();
    }
    for(const [point,fill] of [[pts[0],'#f3d498'],[pts.at(-1),'#b4cde6']]) {
      ctx.fillStyle=fill;ctx.beginPath();ctx.arc(point[0],point[1],4,0,2*Math.PI);ctx.fill();
    }
    const target=Math.max(1e-20,.16*width/scale),base=10**Math.floor(Math.log10(target));
    const length=([1,2,5,10].filter(v=>v*base<=target).at(-1)||1)*base;
    ctx.fillStyle='#e3f1e4';ctx.strokeStyle='#c3daca';ctx.lineWidth=2;
    const sx=14,sy=height-12,dx=length*scale;
    ctx.beginPath();ctx.moveTo(sx,sy-4);ctx.lineTo(sx,sy);ctx.lineTo(sx+dx,sy);ctx.lineTo(sx+dx,sy-4);ctx.stroke();
    ctx.font='13px system-ui';ctx.fillText(F(length,3)+' nm',sx,sy-7);
  }
  function stop(){clearInterval(timer);timer=null;$('play').textContent='Örnekleri oynat';$('play').setAttribute('aria-pressed','false');}
  function drawChain() {
    const key=state.rg+'|'+seed;
    if(coordsKey!==key){xyz=phys.gaussian(state.rg,seed,64,384);coordsKey=key;}
    $('sampleRg').textContent='Rg = '+F(xyz.sampleRg_nm)+' nm';
    $('ensembleRg').textContent='√⟨Rg²⟩ = '+F(xyz.ensembleRmsRg_nm)+' nm';
    $('geometry').textContent='Ölçek '+F(state.rg)+' nm hedef yarıçapa göre belirlenir. Gösterilen zincir 64 Gauss adımıdır. Topluluk ortalaması için ayrıca 384 zincir çekilir. Döndürmek için sürükle veya ok tuşlarını kullan.';
    coilPaint();
  }
  function nextFrame(){seed=(seed+1)>>>0;drawChain();}
  $('next').onclick=()=>{stop();nextFrame();};
  $('play').onclick=()=>{if(timer){stop();return;}$('play').textContent='Oynatmayı durdur';$('play').setAttribute('aria-pressed','true');nextFrame();timer=setInterval(nextFrame,900);};
  const canvas=$('coil');
  canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);};
  canvas.onpointermove=e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*.012;pitch+=(e.clientY-drag[1])*.012;drag=[e.clientX,e.clientY];coilPaint();};
  canvas.onlostpointercapture=()=>{drag=null;};
  canvas.onpointerup=()=>{drag=null;};
  canvas.onkeydown=e=>{const ds={ArrowLeft:[-.12,0],ArrowRight:[.12,0],ArrowUp:[0,-.12],ArrowDown:[0,.12]}[e.key];if(ds){e.preventDefault();yaw+=ds[0];pitch+=ds[1];coilPaint();}};
  $('reset').onclick=()=>{yaw=.3;pitch=-.3;coilPaint();};
  window.addEventListener('resize',coilPaint);
  document.addEventListener('visibilitychange',()=>{if(document.hidden)stop();});
  function ratioPlot(r) {
    const svg=$('regimePlot');svg.replaceChildren();
    const X=x=>32+408*(x+2)/4;
    add(svg,'rect',{x:X(-2),y:48,width:204,height:24,fill:'#daeade',rx:4});
    add(svg,'rect',{x:X(0),y:48,width:204,height:24,fill:'#cedee2',rx:4});
    add(svg,'line',{x1:X(0),x2:X(0),y1:34,y2:83,stroke:'#3e6047','stroke-dasharray':'3 3'});
    for(let i=-2;i<=2;i++)add(svg,'text',{x:X(i),y:95,'text-anchor':'middle',fill:'#415c4c','font-size':14},F(10**i));
    add(svg,'text',{x:236,y:123,'text-anchor':'middle',fill:'#385543','font-size':15},'Derişim oranı · c/c*');
    const pos=Math.min(2,Math.max(-2,Math.log10(r)));
    add(svg,'circle',{cx:X(pos),cy:59,r:8,fill:'#183c2f',stroke:'white','stroke-width':2});
    add(svg,'text',{x:X(pos),y:27,'text-anchor':pos>1.45?'end':pos<-1.45?'start':'middle',fill:'#183c2f','font-size':16,'font-weight':'600'},'c/c* = '+F(r,3));
  }
  function refresh(redrawBench=false) {
    const mode=$('mode').value;
    $('sampleField').hidden=mode!=='measured';$('mwField').hidden=mode==='measured';$('rgField').hidden=mode!=='custom';$('fitField').hidden=mode!=='model';
    const rec=mode==='measured'?all.find(x=>x.id===$('sample').value):null;
    $('context').textContent=mode==='custom'?'Kendi numunen için seyreltik çözeltide, aynı çözücü ve sıcaklıkta ölçülmüş yarıçapı ve numunenin M𝑤 değerini gir.':catalog.polymer+' · '+catalog.solvent+' · '+F(catalog.T_K-273.15,4)+' °C · seyreltik çözeltide θ referansı';
    const mw=rec?.Mw_g_mol??$('mw').valueAsNumber*1000;
    const rg=rec?.Rg_nm??(mode==='model'&&mw>0?phys.predict(mw,$('fit').value==='free'?free:theta):$('rg').valueAsNumber);
    const c=$('c').valueAsNumber;
    const bad=!(mw>0 && Number.isFinite(mw))?'M𝑤 değerini kg/mol cinsinden gir.':!(rg>0 && Number.isFinite(rg))?'Bu numunenin seyreltik çözelti Rg ölçümünü nm cinsinden gir.':!(c>=0 && Number.isFinite(c))?'Polimer derişimi sıfır veya pozitif bir sayı olmalı.':'';
    state=null;result=null;
    $('invalid').hidden=!bad;$('invalid').textContent=bad;$('computed').hidden=!!bad;
    const pair=new URLSearchParams();if(rec&&rows.includes(rec))pair.set('sample',rec.id);
    if(c>=0&&Number.isFinite(c))pair.set('c',String(c));
    $('pairLink').href='chain-solvents.html?'+pair;
    const handoff=new URLSearchParams({source:'chain',system:mode==='custom'?'custom':'ps'});
    if(mode!=='custom'){handoff.set('solvent','cyclohexane');handoff.set('T',String(catalog.T_K));}
    if(rec)handoff.set('sample',rec.id);
    if(Number.isFinite(mw)&&mw>0)handoff.set('Mw',String(mw));
    if(Number.isFinite(c)&&c>=0)handoff.set('c',String(c));
    if(mode==='custom'&&query.has('name'))handoff.set('name',query.get('name'));
    if(query.has('phaseCount'))handoff.set('phaseCount',query.get('phaseCount'));
    $('viscosityLink').href='viscosity.html?'+handoff;
    if(bad){stop();if(redrawBench)updateEvidence();return;}
    const estimate=phys.overlap(mw,rg,c);
    if (![estimate.cstar_g_L,estimate.ratio].every(Number.isFinite)||estimate.cstar_g_L<=0){
      $('computed').hidden=true;$('invalid').hidden=false;$('invalid').textContent='Bu büyüklükler sayısal hesap aralığının dışında.';stop();return;
    }
    state={mw,rg,c,mode,rec};result=estimate;
    $('rgVal').textContent=F(rg)+' nm';
    $('rgSrc').textContent=rec?scaleName(rec)+' · '+rec.id:mode==='custom'?'Kullanıcının deneysel girdisi':(mw<free.minMw||mw>free.maxMw?'Uyum aralığı dışında tahmin':'Uyum aralığı içinde tahmin');
    $('cstar').textContent=F(estimate.cstar_g_L)+' g/L';$('ratio').textContent=F(estimate.ratio);
    const r=estimate.ratio;
    $('regime').textContent=r===0?'Saf çözücü sınırı':r<1?'Seçilen c* değerinin altında':r>1?'Seçilen c* değerinin üzerinde':'Örtüşme ölçeğinde';
    $('ratioExplain').textContent='c = '+F(c)+' g/L için zincirlerin seyreltik çözeltideki eşdeğer küre hacimlerinin toplamı, çözelti hacminin '+F(r)+' katıdır. c/c* büyüdükçe farklı zincirlerin aynı hacmi paylaşması beklenir.';
    const pm=rec?.Rg_reported_pm_nm;
    $('radiusSensitivity').hidden=!pm;
    if(pm){const low=phys.overlap(mw,rg+pm,0).cstar_g_L,high=phys.overlap(mw,rg-pm,0).cstar_g_L;
      $('radiusSensitivity').textContent='Sabit M𝑤 ile kaynakta verilen Rg ± '+F(pm)+' nm, c* için '+F(low)+'–'+F(high)+' g/L verir. Bu aralık yalnız yarıçap duyarlılığıdır.';}
    $('ratioSlide').value=String(Math.min(2,Math.max(-2,Math.log10(r))));ratioPlot(r);
    $('chainPanel').hidden=mode==='custom'||!catalog.thetaReference;
    if($('chainPanel').hidden)stop();else drawChain();
    if(redrawBench)updateEvidence();
  }
  const scopeLabel = x => x==='interpolation'?'Uyum aralığında':'Uyum aralığı dışında';
  function massPlot(group) {
    const isExt=group==='external', shown=isExt?all:rows, axis=$('massPlot');axis.replaceChildren();
    const lo= Math.min(...shown.map(x=>x.Mw_g_mol))/1.4, hi=Math.max(...shown.map(x=>x.Mw_g_mol))*1.25;
    const lowY=Math.min(...shown.map(x=>x.Rg_nm),phys.predict(lo,free),phys.predict(lo,theta))*.80;
    const highY=Math.max(...shown.map(x=>x.Rg_nm),phys.predict(hi,free),phys.predict(hi,theta))*1.18;
    const left=56,right=548,top=26,bottom=245;
    const X=x=>left+(right-left)*Math.log(x/lo)/Math.log(hi/lo);
    const Y=y=>bottom-(bottom-top)*Math.log(y/lowY)/Math.log(highY/lowY);
    const text=(x,y,label,attrs={})=>add(axis,'text',{x,y,'font-size':15,fill:'#3e594a',...attrs},label);
    add(axis,'rect',{x:X(free.minMw),y:top,width:X(free.maxMw)-X(free.minMw),height:bottom-top,fill:'#eff5ee'});
    for(let y of logAxisTicks(lowY,highY,5)){
      add(axis,'line',{x1:left,x2:right,y1:Y(y),y2:Y(y),stroke:'#dae2da'});text(left-7,Y(y)+4,F(y,3),{'text-anchor':'end'});
    }
    for(let x of logAxisTicks(lo/1000,hi/1000,5)){
      const p=X(x*1000);add(axis,'line',{x1:p,x2:p,y1:bottom,y2:bottom+5,stroke:'#8da194'});
      text(p,bottom+23,F(x,3),{'text-anchor':'middle'});
    }
    add(axis,'path',{d:'M'+left+' '+top+'V'+bottom+'H'+right,stroke:'#8da194',fill:'none'});
    text(left,16,'Rg (nm)');text((left+right)/2,298,'M𝑤 (kg/mol) · logaritmik eksenler',{'text-anchor':'middle'});
    for(const [fit,color,dash] of [[free,'#1b7459','none'],[theta,'#795b78','6 5']]){
      add(axis,'line',{x1:X(lo),x2:X(hi),y1:Y(phys.predict(lo,fit)),y2:Y(phys.predict(hi,fit)),stroke:color,'stroke-dasharray':dash,'stroke-width':2.4});
    }
    for(const row of shown){
      const x=X(row.Mw_g_mol),y=Y(row.Rg_nm),a=ext.includes(row);
      if(row.id===state?.rec?.id)add(axis,'circle',{cx:x,cy:y,r:9,fill:'none',stroke:'#98732c','stroke-width':1.7});
      if(row.Rg_reported_pm_nm) {
        const e=row.Rg_reported_pm_nm;add(axis,'line',{x1:x,x2:x,y1:Y(row.Rg_nm-e),y2:Y(row.Rg_nm+e),stroke:'#253f32','stroke-width':1.5});
      }
      const mark=a?add(axis,'path',{d:'M'+x+' '+(y-5)+'l5 5 -5 5 -5 -5z',fill:'white',stroke:'#94582d','stroke-width':2}):add(axis,'circle',{cx:x,cy:y,r:4.1,fill:'#233d2e'});
      const title=document.createElementNS(ns,'title');title.textContent=(row.label||row.id)+' · '+F(row.Rg_nm)+' nm';mark.append(title);
    }
    $('externalKey').hidden=!isExt;
  }
  function updateEvidence() {
    const group=$('validation').value;
    const isExt=group==='external'&&external, b=isExt?external:loo;
    $('evidenceSystem').textContent='Polistiren / siklohekzan · 34,5 °C';
    $('benchDescription').textContent=isExt?'Altı yüksek mol kütleli fraksiyona uydurulan iki model ayrı yayındaki iki düşük kütleli SAXS ölçümüne uygulanır.':'Her turda bir fraksiyon ayrılır; iki modelin katsayıları kalan fraksiyonlardan yeniden hesaplanır. Hatalar ayrılan ölçümlerde hesaplanır.';
    $('errFree').textContent='%'+F(b.metrics.free.rmse_pct,3);
    $('errFixed').textContent='%'+F(b.metrics.theta.rmse_pct,3);
    $('nu').textContent=F(free.nu,4);
    $('heldNote').textContent='Hata, her fraksiyonun Rg değerine göre hesaplanan bağıl hataların kareli ortalamasının kareköküdür. Yeşil bant tüm uyum serisinin M𝑤 aralığını, çizgiler tüm uyum serisinden belirlenen modelleri gösterir. '+(isExt?'Dış ölçümler hem yöntem hem mol kütlesi bakımından farklıdır. İki nokta genel bir doğruluk oranı belirlemeye yetmez.':'Tablo her fraksiyonun kendi dışarıda bırakma tahminini verir. Bu altı fraksiyon aynı yayındaki iki ana örnekten gelir; komşu kütleler bağımsız yayınlardan alınmamıştır.');
    const body=$('records');body.replaceChildren();
    for(const r of b.folds){
      const row=all.find(q=>q.id===r.id),tr=body.insertRow();
      const button=document.createElement('button');button.type='button';button.className='rowselect';button.textContent=row.label||row.id;
      button.setAttribute('aria-label',(row.label||row.id)+' fraksiyonunu seç');
      button.onclick=()=>{stop();$('mode').value='measured';$('sample').value=row.id;refresh(true);$('sample').focus();};
      tr.insertCell().append(button);
      for(const val of [F(r.Mw_g_mol/1000),F(r.observed_nm),F(r.free_nm),F(r.theta_nm),scopeLabel(r.domain)])tr.insertCell().textContent=val;
    }
    massPlot(isExt?'external':'loo');
  }
  $('downloadRows').onclick=()=>{
    const key=$('validation').value, b=key==='external'&&external?external:loo;
    const src=key==='external'?other.source:catalog.source;
    exportCsv('PS_siklohekzan_Rg_'+key+'.csv',[
      ['sample','source_doi','split','T_K','Mw_g_mol','observed_Rg_nm','predicted_free_Rg_nm','predicted_theta_Rg_nm','domain','training_ids'],
      ...b.folds.map(r=>[r.id,src.doi,key,catalog.T_K,r.Mw_g_mol,r.observed_nm,r.free_nm,r.theta_nm,r.domain,r.train_ids.join('|')])
    ]);
  };
  $('downloadResult').onclick=()=>{
    if(!state||!result)return;
    const rec=state.rec,src=state.mode==='custom'?null:rec&&ext.includes(rec)?other.source:catalog.source;
    const temp = state.mode==='custom'?(T??''):catalog.T_K;
    exportCsv('Zincir_boyutu_ortusme.csv',[
      ['sample','basis','solvent','temperature_K','reference_doi','Mw_g_mol','Rg_nm','c_g_L','cstar_g_L','c_over_cstar','definition','model_nu'],
      [rec?.id||(state.mode==='custom'?'user':'PS-theta-prediction'),state.mode,state.mode==='custom'?'':catalog.solvent,temp,src?.doi??'',state.mw,state.rg,state.c,result.cstar_g_L,result.ratio,'3Mw/(4π NA Rg^3)',state.mode==='model'?($('fit').value==='free'?free.nu:theta.nu):'']
    ]);
  };
  $('mode').onchange=()=>{
    stop();if($('mode').value==='model'&&!$('mw').value)$('mw').value=String((all.find(x=>x.id===$('sample').value)||rows[0]).Mw_g_mol/1000);
    refresh(true);
  };
  $('sample').onchange=()=>{stop();refresh(true);};
  $('mw').oninput=()=>{stop();refresh();};$('rg').oninput=()=>{stop();refresh();};
  $('fit').onchange=()=>{stop();refresh();};$('c').oninput=()=>refresh();
  $('ratioSlide').oninput=()=>{if(result){$('c').value=Number((10**Number($('ratioSlide').value)*result.cstar_g_L).toPrecision(12));refresh();}};
  $('validation').onchange=updateEvidence;
  $('methodText').textContent='Referans seri 34,5 °C’de siklohekzanda ölçülen altı lineer PS fraksiyonudur. Kuvvet yasaları ln Rg ve ln M𝑤 uzayında her fraksiyona eşit ağırlık verilerek uydurulur. Rg, bir zincirin gözlenen anlık geometrisi değil, saçılmada elde edilen ortalama kare boyutun kareköküdür.';
  const refs=[catalog.source,...(other?[other.source]:[])];
  if(catalog.source.numericCheck)refs.push({label:catalog.source.numericCheckLocation+' · aynı ölçümlerin sayısal çapraz kontrolü',url:catalog.source.numericCheck});
  for (const s of refs){const li=document.createElement('li'),a=document.createElement('a');a.textContent=(s.label||s.caption||s.doi);a.href=s.url;li.append(a);$('sources').append(li);}
  $('app').hidden=false;
  refresh(true);
})();
