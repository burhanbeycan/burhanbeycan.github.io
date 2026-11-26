(() => {
  'use strict';
  const el=id=>document.getElementById(id);
  const theta=window.THETA_FAMILY_BENCHMARK, hydro=window.HYDRO_BENCHMARK;
  const tr=(value,digits=4)=>Number.isFinite(value)?value.toLocaleString('tr-TR',{maximumSignificantDigits:digits}):'—';
  const colors=['#196452','#9b532f','#426a98','#8a4e7c','#79702a','#376b6c'];
  const names={
    pbd:'Polibütadien',pb:'Polibütadien',pi:'Poliizopren',pib:'Poliizobütilen',ps:'Polistiren',
    pams:'Poli(α-metilstiren)',palphaMS:'Poli(α-metilstiren)',diene:'Polidienler',polydiene:'Polidienler',
    polydienes:'Polidienler',saturated:'Doymuş omurgalı polimerler',polyolefin:'Doymuş poliolefinler',
    vinylaromatic:'Vinil aromatikler','ps-benzene':'Polistiren · benzen','pullulan-water':'Pullulan · su',
    cyclohexane:'Siklohekzan',benzene:'Benzen',water:'Su',pullulan:'Pullulan','water-NaN3-0.2wt':'Su (%0,2 kütlece NaN₃)', 'isoamyl isovalerate':'İzoamil izovalerat','1,4-dioxane':'1,4-Dioksan'
  };
  const label=x=>names[x]||names[String(x).toLowerCase()]||String(x||'—');
  const numeric=x=>x===null||x===undefined||x===''?null:Number.isFinite(+x)?+x:null;
  const option=(select,value,text)=>select.add(new Option(text,value));
  const ns='http://www.w3.org/2000/svg';
  const svgItem=(node,name,attrs,content)=>{
    const tag=document.createElementNS(ns,name);
    for(const [key,value] of Object.entries(attrs))tag.setAttribute(key,String(value));
    if(content!=null)tag.textContent=content;node.append(tag);return tag;
  };
  const link=(parent,url,text)=>{
    if(!url)return;
    const a=document.createElement('a');a.href=url;a.textContent=text;parent.append(a);
  };
  const methods={
    theta:['Teta · mol kütlesi','Teta · yapı tanımlayıcılı Ridge'],
    hydro:['Mol kütlesi · Mark–Houwink','Ölçülen Rɡ · Flory–Fox']
  };
  const protocolNames={chemistry:'Bir polimer kimyasını ayır',backbone:'Bir omurga ailesini ayır',external:'Farklı deney koşuluna aktar'};
  function packTheta(){
    if(!theta?.families?.length)return [];
    return theta.families.map(f=>({
      key:f.family,title:f.label||label(f.family),split:'backbone',nTrain:f.nTrain,
      train:f.trainingFamilies||[],range:[numeric(f.trainMassMin),numeric(f.trainMassMax)],
      rows:f.folds.map(r=>({
        id:r.id,polymer:r.chemistry,group:f.family,observed:r.observed,a:r.theta,b:r.ridge,
        mass:r.mass,rg:null,solvent:r.solvent,temp:numeric(r.T_C),domain:r.domain,
        descDomain:r.descriptorDomain,ref:r.ref,doi:r.source_doi,table:r.table,source:r.source_url||null,
        basis:r.mass_basis||null,tempBasis:r.temperature_basis
      }))
    }));
  }
  function packHydro(){
    if(!hydro?.folds?.length)return [];
    return hydro.folds.map(f=>({
      key:f.group,title:f.label||label(f.group),split:f.split,nTrain:f.training_n,
      train:f.trainFamilies||[],range:[numeric(f.training_interval?.M_min),numeric(f.training_interval?.M_max)],u:f.u_fit,
      rows:f.rows.map(r=>({
        id:r.id,polymer:r.chemistry||r.family,group:f.group,observed:r.eta_mL_g,a:r.mass_prediction,b:r.size_prediction,
        mass:r.M_g_mol,rg:r.Rg_nm,solvent:r.solvent,temp:Number.isFinite(r.T_K)?r.T_K-273.15:numeric(r.T_C),
        domain:r.mass_domain,rgDomain:r.radius_domain,ref:r.original_ref,doi:r.source_doi,rgDoi:r.radius_doi,table:r.table,source:r.source_url||null,basis:r.mass_basis||null
      }))
    }));
  }
  const groups={theta:packTheta(),hydro:packHydro()};
  const q=new URLSearchParams(location.search);
  let cohort=q.get('cohort')==='hydro'&&groups.hydro.length?'hydro':groups.theta.length?'theta':'hydro';
  let activeGroups=[],activeRows=[],record=null;
  const fail=msg=>{el('familyError').textContent=msg;el('familyError').hidden=false;};
  if(!groups.theta.length&&!groups.hydro.length){fail('Sınama sonuçları yüklenemedi. Çevrimdışı paket kullanıyorsan tüm dosyaları aynı dizine çıkar.');return;}
  if(!groups.theta.length||!groups.hydro.length)fail('Bir sınamanın sonuç dosyası yüklenemedi. Şu anda erişilebilen deneyi inceleyebilirsin.');
  el('tabTheta').disabled=!groups.theta.length;el('tabHydro').disabled=!groups.hydro.length;
  const selectedMethod=()=>el('familyMethod').value==='b'?'b':'a';
  const available=()=>groups[cohort].filter(g=>cohort==='theta'||g.split===el('familyProtocol').value);
  function updateTabs(initial=false){
    const isTheta=cohort==='theta';
    el('tabTheta').setAttribute('aria-pressed',String(isTheta));el('tabHydro').setAttribute('aria-pressed',String(!isTheta));
    const protocol=el('familyProtocol');
    protocol.replaceChildren();
    if(!isTheta){
      const splits=new Set(groups.hydro.map(g=>g.split));
      for(const [key,name] of Object.entries(protocolNames))if(splits.has(key))option(protocol,key,name);
      if(splits.has(initial?q.get('split'):'backbone'))protocol.value=initial?q.get('split'):'backbone';
      else if(splits.has('backbone'))protocol.value='backbone';
    }
    el('protocolField').hidden=isTheta||protocol.options.length<2;
    const method=el('familyMethod');method.replaceChildren();
    methods[cohort].forEach((name,i)=>option(method,i?'b':'a',name));method.value=initial&&q.get('method')==='a'?'a':'b';
    updateGroupOptions(initial);
  }
  function updateGroupOptions(initial=false){
    const list=available(),fold=el('familyFold');
    fold.replaceChildren();if(list.length>1)option(fold,'all','Tüm saklı gruplar');
    for(const g of list)option(fold,g.key,g.title);
    if(initial&&list.some(g=>g.key===q.get('family')))fold.value=q.get('family');
    updateGroup();
  }
  function updateGroup(){
    const all=available(),choice=el('familyFold').value;
    activeGroups=choice==='all'?all:all.filter(g=>g.key===choice);
    activeRows=activeGroups.flatMap(g=>g.rows.map(r=>({...r,fold:g})));
    const rec=el('familyRecord');rec.replaceChildren();
    activeRows.forEach((r,i)=>option(rec,String(i),label(r.polymer)+' · '+tr(r.mass/1000)+' kg/mol · '+r.id));
    draw();
  }
  function mape(rows,key){
    return rows.reduce((n,r)=>n+Math.abs(r[key]/r.observed-1),0)*100/rows.length;
  }
  function updateText(){
    const list=available(),allCount=list.reduce((v,g)=>v+g.rows.length,0);
    el('familyScope').textContent=activeGroups.length===1?activeGroups[0].title:activeGroups.length+' saklı grup · '+activeRows.length+' ölçüm';
    const isTheta=cohort==='theta', protocol=el('familyProtocol').value;
    el('familySensitivity').hidden=!isTheta;
    el('familyStudy').textContent=isTheta
      ? 'Teta çözücülerindeki fraksiyonlar için iki model de aynı mol kütlesini kullanır. İkinci model ayrıca tekrar biriminin kimyasal yapısından çıkarılan tanımlayıcıları kullanır. Çözücü ve sıcaklık polimere göre değişir.'
      : protocol==='external'
        ? 'Referans polimerlerde öğrenilen katsayılar başka bir çözücüdeki ölçümlere uygulanır. Bu aktarımda polimer ile deney koşulu birlikte değişebilir.'
        : 'Siklohekzandaki aynı fraksiyonlar üzerinde mol kütlesi tahminini, hedef numunenin saçılmayla ölçülmüş jirasyon yarıçapını da kullanan modelle karşılaştır. Derlemenin referans sıcaklığı 25 °C’dir. Bazı satırlarda yalnız bu nominal sıcaklık bildirilmiştir.';
    el('familyPhysics').textContent=isTheta
      ? 'Teta koşulu için her iki model [η] ∝ M¹ᐟ² bağıntısını kullanır. Birincisi tüm eğitim polimerlerinden ortak bir ön çarpan öğrenir. Ridge modeli bu ön çarpanı, tekrar birimi başına omurga bağ sayısının tekrar birimi mol kütlesine oranının logaritmasına bağlar. Tanımlayıcı eğitimde standartlaştırılır ve ceza katsayısı λ = 1 olarak sabitlenir. Her turda katsayılar yalnız eğitim gruplarından uyarlanır.'
      : 'Mark–Houwink modeli log [η] ile log M arasındaki ilişkiyi öğrenir. Flory–Fox modelinde ölçülen Rɡ³/M boyut ölçeği kullanılır. Katsayı eğitim polimerlerinden alınır. Saklı numunenin [η] ölçümü her iki modelde de yalnız değerlendirmede kullanılır.';
    el('familyLimits').textContent=isTheta
      ? 'Teta koşulunda farklı zincirler farklı segment uzunluklarına ve yerel konformasyonlara sahip olabilir. Bu küçük aile kümesi, yeni bir kimya için genel doğruluk veya güven aralığı belirlemez. Mol kütlesinin eğitim aralığına girmesi, kimyasal benzerliği garanti etmez.'
      : 'Zincir boyutlu tahmin hedef numunede ayrıca deneysel Rɡ gerektirir. Saçılmanın ve viskozimetrinin hassas olduğu dağılım ortalamaları farklıdır. Sınırlı kimya sayısı ile farklı çözücü ve laboratuvarlara aktarım ayrı sınamalar gerektirir.';
    el('familyErrorIntro').textContent=list.length+' grup ve '+allCount+' saklı ölçüm. Grupların eşit ağırlıklı ortalaması az ölçümlü aileleri görünür kılar. Birleşik ortalamada her ölçüm eşit ağırlık alır.';
    const g=activeGroups.length===1?activeGroups[0]:null;
    const training=g&&Array.isArray(g.train)&&g.train.length?'Eğitim grupları · '+g.train.map(label).join(', ')+'. ':'';
    const range=g&&g.range.every(Number.isFinite)?'Eğitim kütle aralığı '+tr(g.range[0]/1000)+'–'+tr(g.range[1]/1000)+' kg/mol.':'';
    const ood=activeRows.filter(r=>r.descDomain==='extrapolation').length;
    const supportNote=isTheta?' Yapı tanımlayıcısı '+activeRows.length+' gözlemin '+ood+' tanesinde ilgili eğitim aralığının dışındadır.':'';
    el('familyFitInfo').textContent=(g?training+(g.nTrain?g.nTrain+' eğitim ölçümü. ':'')+range:'Her işaret, kendi grubunun dışarıda tutulduğu uyumdan gelir. Sonuçları grup seçerek inceleyebilirsin.')+supportNote;
  }
  function updateScores(){
    const all=available(), tbody=el('familyScores');tbody.replaceChildren();
    const labels=methods[cohort];el('familyMethodA').textContent=labels[0]+' · MAPE';el('familyMethodB').textContent=labels[1]+' · MAPE';
    const cells=(title,n,x,y,style)=>{
      const row=tbody.insertRow();row.className=style;
      const th=document.createElement('th');th.scope='row';th.textContent=title;row.append(th);
      [n,'%'+tr(x,3),'%'+tr(y,3)].forEach(s=>row.insertCell().textContent=s);
    };
    for(const g of all)cells(g.title,g.rows.length,mape(g.rows,'a'),mape(g.rows,'b'),el('familyFold').value===g.key?'family-focus':'');
    if(all.length>1){
      const rows=all.flatMap(g=>g.rows);
      const macro=k=>all.reduce((v,g)=>v+mape(g.rows,k),0)/all.length;
      cells('Gruplar eşit ağırlıklı',all.length+' grup',macro('a'),macro('b'),'family-summary');
      cells('Ölçümler eşit ağırlıklı',rows.length,mape(rows,'a'),mape(rows,'b'),'family-summary');
    }
  }
  function detail(){
    const key=selectedMethod();record=activeRows[+el('familyRecord').value];
    if(!record)return;
    const r=record,bias=100*(r[key]/r.observed-1);
    el('familyMass').textContent=tr(r.mass/1000)+' kg/mol';
    el('familyObs').textContent=tr(r.observed)+' mL/g';el('familyPred').textContent=tr(r[key])+' mL/g';
    el('familyRgRow').hidden=!Number.isFinite(r.rg);el('familyRg').textContent=tr(r.rg)+' nm';
    el('familyBias').textContent=(bias>0?'+':'')+tr(bias,3)+' %';
    el('familyCondition').textContent=[label(r.polymer),r.solvent?label(r.solvent):'',Number.isFinite(r.temp)?(r.tempBasis==='approximate_compilation'?'≈ ':'')+tr(r.temp)+' °C'+(cohort==='hydro'&&r.solvent==='cyclohexane'?' (nominal)':''):''].filter(Boolean).join(' · ');
    const massDom=['inside','interpolation','within'].includes(r.domain)?'Mol kütlesi eğitim aralığında.':['outside','extrapolation','below','above'].includes(r.domain)?'Mol kütlesi eğitim aralığının dışında.':'';
    const basis={M_abs:'derlemedeki mutlak mol kütlesi',Mw:'ağırlık ortalamalı M𝑤',M_table:'derlemede verilen mol kütlesi'};
    const addDomain=r.descDomain==='extrapolation'?'Yapı tanımlayıcısı eğitim aralığının dışında.':r.descDomain==='interpolation'?'Yapı tanımlayıcısı eğitim aralığında.':r.rgDomain==='outside'?'Rɡ eğitim aralığının dışında.':r.rgDomain==='inside'?'Rɡ eğitim aralığında.':'';
    el('familyDomain').textContent=[massDom,addDomain,r.basis?'Kütle tanımı · '+(basis[r.basis]||r.basis)+'.':''].filter(Boolean).join(' ');
    const src=el('familyReference');src.replaceChildren();
    const refs=r.ref?'Özgün yayın · derleme kaynağı '+r.ref:r.table?'Viskozite ölçümü · Tablo '+r.table:'Ölçüm kaynağı';
    const url=r.doi?(String(r.doi).startsWith('https://')?r.doi:'https://doi.org/'+r.doi):r.source;
    if(url)link(src,url,refs||'Ölçüm kaynağı');else if(refs)src.textContent=refs;else src.textContent='';
    if(r.rgDoi&&r.rgDoi!==r.doi){src.append(' · ');link(src,'https://doi.org/'+r.rgDoi,'Saçılma kaynağı');}
    for(const button of el('familyRows').querySelectorAll('button[data-index]'))button.closest('tr').classList.toggle('family-focus',+button.dataset.index===+el('familyRecord').value);
  }
  function symbol(svg,x,y,style,color,r){
    const props={fill:'#ffffff','fill-opacity':.83,stroke:color,'stroke-width':1.85};
    switch(style%4){
      case 0:return svgItem(svg,'circle',{...props,cx:x,cy:y,r});
      case 1:return svgItem(svg,'path',{...props,d:'M'+x+' '+(y-r-1)+'l'+(r+1)+' '+(r+1)+' '+(-r-1)+' '+(r+1)+' '+(-r-1)+' '+(-r-1)+'z'});
      case 2:return svgItem(svg,'rect',{...props,x:x-r,y:y-r,width:2*r,height:2*r});
      default:return svgItem(svg,'path',{...props,d:'M'+x+' '+(y-r-1)+'l'+(r+1)+' '+(2*r+1)+'h'+(-2*r-2)+'z'});
    }
  }
  function drawChart(){
    const svg=el('familyParity');svg.replaceChildren();
    const key=selectedMethod(), values=activeRows.flatMap(r=>[r.observed,r[key]]).filter(v=>Number.isFinite(v)&&v>0);
    if(!values.length)return;
    const min=Math.min(...values),max=Math.max(...values),pad=Math.max(.08,(Math.log10(max)-Math.log10(min))*.065);
    const lo=Math.log10(min)-pad,hi=Math.log10(max)+pad, left=94,top=20,w=405,h=261;
    const X=v=>left+w*(Math.log10(v)-lo)/(hi-lo), Y=v=>top+h-h*(Math.log10(v)-lo)/(hi-lo);
    const tick=(x,y,s,align='middle')=>svgItem(svg,'text',{x,y,'text-anchor':align,'font-size':14},s);
    const tickCandidates=[];
    for(let d=Math.floor(lo);d<=Math.ceil(hi);d++)for(const v of [1,2,5]){const n=v*10**d;if(Math.log10(n)>=lo&&Math.log10(n)<=hi)tickCandidates.push(n);}
    let ticks=tickCandidates.filter((_,i)=>i%Math.ceil(tickCandidates.length/6)===0);
    if(ticks.length<2)ticks=[min,max];
    for(const v of [...new Set(ticks)]){
      const x=X(v),y=Y(v);
      svgItem(svg,'path',{d:'M'+x+' '+top+'v'+h+' M'+left+' '+y+'h'+w,fill:'none',stroke:'#e0e8df','stroke-width':1});
      tick(x,top+h+21,tr(v));tick(left-9,y+4,tr(v),'end');
    }
    svgItem(svg,'path',{d:'M'+left+' '+top+'v'+h+'h'+w,stroke:'#84998c',fill:'none'});
    svgItem(svg,'line',{x1:left,y1:top+h,x2:left+w,y2:top,stroke:'#7b8981','stroke-width':1.6,'stroke-dasharray':'5 5'});
    tick(left+w/2,top+h+55,'Ölçülen [η] (mL/g)');
    svgItem(svg,'text',{transform:'translate(28 '+(top+h/2)+') rotate(-90)','text-anchor':'middle','font-size':16},'Tahmin edilen [η] (mL/g)');
    const all=available(),legend=el('familyLegend');legend.replaceChildren();
    for(const g of activeGroups){
      const i=all.indexOf(g),c=colors[i%colors.length],span=document.createElement('span'),mark=document.createElementNS(ns,'svg');
      mark.setAttribute('viewBox','0 0 14 14');mark.setAttribute('aria-hidden','true');mark.style.cssText='width:14px;height:14px;min-height:0;background:none';symbol(mark,7,7,i,c,3.5);span.append(mark,document.createTextNode(g.title));legend.append(span);
    }
    activeRows.forEach((r,i)=>{
      const j=all.indexOf(r.fold),x=X(r.observed),y=Y(r[key]),dot=symbol(svg,x,y,j,colors[j%colors.length],3.8);
      svgItem(dot,'title',{},label(r.polymer)+' · '+r.id+' · ölçüm '+tr(r.observed)+' · model '+tr(r[key])+' mL/g');
      dot.style.cursor='pointer';dot.addEventListener('click',()=>selectRecord(i));
    });
    if(record)svgItem(svg,'circle',{cx:X(record.observed),cy:Y(record[key]),r:8,class:'active','pointer-events':'none'});
  }
  function drawRows(){
    const body=el('familyRows');body.replaceChildren();
    el('detailA').textContent=methods[cohort][0]+' [η]';el('detailB').textContent=methods[cohort][1]+' [η]';
    activeRows.forEach((r,i)=>{
      const trEl=body.insertRow(),button=document.createElement('button');
      button.type='button';button.dataset.index=i;button.textContent=r.id;button.addEventListener('click',()=>selectRecord(i));
      trEl.insertCell().append(button);
      for(const s of [label(r.polymer),tr(r.mass/1000),tr(r.observed),tr(r.a),tr(r.b)])trEl.insertCell().textContent=s;
    });
  }
  function selectRecord(i){el('familyRecord').value=String(i);detail();drawChart();}
  function draw(){updateText();updateScores();drawRows();detail();drawChart();}
  function csv(){
    const headings=['cohort','split','held_group','measurement_id','chemistry','source_doi','source_ref','temperature_C','solvent','M_g_mol','mass_basis','Rg_nm','eta_observed_mL_g','baseline_mL_g','second_model_mL_g','baseline','second_model','mass_domain','descriptor_domain','Rg_domain','Rg_source_doi','training_groups'];
    const values=activeRows.map(r=>[cohort,r.fold.split,r.group,r.id,r.polymer,r.doi,r.ref,r.temp,r.solvent,r.mass,r.basis,r.rg,r.observed,r.a,r.b,...methods[cohort],r.domain,r.descDomain,r.rgDomain,r.rgDoi,r.fold.train.join('|')]);
    const cell=v=>'"'+String(v??'').replaceAll('"','""')+'"';
    const blob=new Blob(['\uFEFF'+[headings,...values].map(row=>row.map(cell).join(',')).join('\r\n')+'\r\n'],{type:'text/csv;charset=utf-8'});
    const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='polimer_aile_sinamasi_'+cohort+'_'+el('familyFold').value+'.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  el('tabTheta').addEventListener('click',()=>{cohort='theta';updateTabs();});
  el('tabHydro').addEventListener('click',()=>{cohort='hydro';updateTabs();});
  el('familyProtocol').addEventListener('change',()=>updateGroupOptions());
  el('familyFold').addEventListener('change',updateGroup);
  el('familyMethod').addEventListener('change',()=>{detail();drawChart();});
  el('familyRecord').addEventListener('change',()=>{detail();drawChart();});
  el('familyCsv').addEventListener('click',csv);
  el('familyStudio').hidden=false;updateTabs(true);
})();
