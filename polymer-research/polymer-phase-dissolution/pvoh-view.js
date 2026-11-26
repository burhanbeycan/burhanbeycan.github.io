(() => {
  'use strict';
  const el=id=>document.getElementById(id);
  const physics=window.PvohViscosity, catalog=window.FAMILY_PVOH;
  if(!physics||!catalog?.rows?.length){el('loadError').hidden=false;return;}
  const {rows}=catalog;
  const result=physics.evaluate(rows);
  const fmt=(v,n=4)=>Number.isFinite(v)?v.toLocaleString('tr-TR',{maximumSignificantDigits:n}):'—';
  const titles={water:'Su',dmso:'DMSO',all:'Su ve DMSO'};
  const color={water:'#20695c',dmso:'#9b552e'};
  const query=new URLSearchParams(location.search);
  const ids=[...new Set(rows.map(r=>r.grade))];
  for(const id of ids)el('grade').add(new Option(id,id));
  el('grade').value=ids.includes(query.get('grade'))?query.get('grade'):'49k-75';
  if(['water','dmso'].includes(query.get('solvent')))el('scope').value=query.get('solvent');
  const incomingC=query.get('c');if(incomingC!==null&&incomingC.trim()&&Number.isFinite(+incomingC)&&+incomingC>=0)el('concentration').value=String(+incomingC);
  function observation(){
    const grade=el('grade').value,c=el('concentration').valueAsNumber;
    const selected=rows.filter(r=>r.grade===grade),g=selected[0];
    el('chemistry').textContent='M𝑤 = '+fmt(g.Mw_g_mol/1000)+' kg/mol · Mₙ = '+fmt(g.Mn_g_mol/1000)+' kg/mol · Đ = '+fmt(g.dispersity)+' · DH ≈ %'+fmt(g.DH_pct)+' · kalan asetat ≈ molce %'+fmt(g.residual_acetate_mol_pct)+'.';
    for(const r of selected){
      el(r.solvent+'Eta').textContent=fmt(r.eta_mL_g)+' mL/g';
      el(r.solvent+'H').textContent='Huggins katsayısı kₕ = '+fmt(r.kH)+'.';
      let note='Derişim için sıfır veya pozitif bir değer gir.';
      if(Number.isFinite(c)&&c>=0){
        const s=physics.relative(r,c);
        note=s.etaRel===null?'[η]c = '+fmt(s.u)+'. Seyreltik hesap aralığının üzerinde.':
          'ηᵣ = '+fmt(s.etaRel,5)+' · [η]c = '+fmt(s.u)+'.';
      }
      el(r.solvent+'Rel').textContent=note;
    }
  }
  function draw(points){
    const svg=el('parity');svg.replaceChildren();
    const add=(name,attrs,txt)=>{
      const item=document.createElementNS('http://www.w3.org/2000/svg',name);
      for(const [k,v] of Object.entries(attrs))item.setAttribute(k,String(v));
      if(txt!=null)item.textContent=txt;svg.append(item);return item;
    };
    const tick=el('scope').value==='water'?25:50;
    const limit=tick*Math.ceil(Math.max(...points.flatMap(r=>[r.observed,r.mass,r.dh]))/tick);
    const X=v=>91+346*v/limit,Y=v=>297-244*v/limit;
    for(let i=0;i<=limit;i+=tick){
      add('line',{x1:X(0),y1:Y(i),x2:X(limit),y2:Y(i),stroke:'#dce4de'});
      add('text',{x:X(0)-9,y:Y(i)+4,'text-anchor':'end',fill:'#405a4b','font-size':15},fmt(i));
      add('text',{x:X(i),y:320,'text-anchor':'middle',fill:'#405a4b','font-size':15},fmt(i));
    }
    add('path',{d:'M'+X(0)+' '+Y(limit)+'V'+Y(0)+'H'+X(limit),stroke:'#779386',fill:'none'});
    add('path',{d:'M'+X(0)+' '+Y(0)+'L'+X(limit)+' '+Y(limit),stroke:'#87958b','stroke-dasharray':'5 4',fill:'none'});
    add('text',{x:X(0),y:23,'font-size':16,fill:'#2e4a3d'},'Tahmin edilen [η] (mL/g)');
    add('text',{x:264,y:351,'text-anchor':'middle','font-size':16,fill:'#2e4a3d'},'Ölçülen [η] (mL/g)');
    for(const r of points){
      const x=X(r.observed),y1=Y(r.mass),y2=Y(r.dh),highlight=r.grade===el('grade').value;
      add('line',{x1:x,x2:x,y1,y2,stroke:color[r.solvent],opacity:.4,'stroke-width':1.3});
      const a=add('circle',{cx:x,cy:y1,r:highlight?6:4.5,fill:'white',stroke:color[r.solvent],'stroke-width':highlight?2.4:1.7});
      const size=highlight?6:4.7;
      const b=add('path',{d:'M'+x+' '+(y2-size)+'l'+size+' '+size+' -'+size+' '+size+' -'+size+' -'+size+'z',fill:color[r.solvent],stroke:'white'});
      for(const [mark,n,p] of [[a,'M𝑤',r.mass],[b,'M𝑤 + DH',r.dh]]){
        const title=document.createElementNS('http://www.w3.org/2000/svg','title');title.textContent=r.grade+' · '+titles[r.solvent]+' · '+n+' · '+fmt(p)+' mL/g';mark.append(title);
      }
    }
  }
  function validation(){
    const scope=el('scope').value, grade=el('grade').value, scores=result.scores[scope];
    const points=result.folds.filter(f=>scope==='all'||f.solvent===scope);
    el('massError').textContent=fmt(scores.mass.rmse)+' mL/g';el('dhError').textContent=fmt(scores.dh.rmse)+' mL/g';
    el('subset').textContent='Beş ayrı numune · '+points.length+' saklı çözücü–numune ölçümü. Hata, saklı noktalardaki farkların kare ortalamasının kareköküdür.';
    el('interpret').textContent='Bağıl RMSE mol kütlesi modelinde %'+fmt(scores.mass.relativeRmse,3)+', hidroliz terimli modelde %'+fmt(scores.dh.relativeRmse,3)+'.';
    const body=el('folds');body.replaceChildren();
    for(const f of points){
      const tr=body.insertRow();if(f.grade===grade)tr.className='highlight';
      const labels=[f.grade,titles[f.solvent],fmt(f.observed),fmt(f.mass),fmt(f.dh),
        f.DH_pct>=f.DH_range[0]&&f.DH_pct<=f.DH_range[1]?'Aralık içinde':'Aralık dışında'];
      for(const l of labels)tr.insertCell().textContent=l;
    }
    const f=points.find(r=>r.grade===grade);
    let msg=grade+' saklandığında dört eğitim numunesi eğime katkı verir. Sabit M𝑤 ve çözücüde +10 yüzde puan DH için [η] çarpanı '+fmt(f.dh10Factor)+'.';
    if(f.Mw_g_mol<f.Mw_range[0]||f.Mw_g_mol>f.Mw_range[1])msg+=' Bu numunenin M𝑤 değeri bu turun eğitim aralığı dışında.';
    el('focusSummary').textContent=msg;
    draw(points);
  }
  el('export').onclick=()=>{
    const header='specimen,solvent,T_K,Mw_g_mol,Mn_g_mol,dispersity,DH_approx_pct,intrinsic_observed_mL_g,kH_from_Huggins,kH_from_Kraemer,held_grade_mass_prediction_mL_g,held_grade_mass_DH_prediction_mL_g,fold_DH10_factor,training_grades,source_doi';
    const lines=result.folds.map(f=>{
      const r=rows.find(v=>v.grade===f.grade&&v.solvent===f.solvent);
      return [r.grade,r.solvent,r.T_K,r.Mw_g_mol,r.Mn_g_mol,r.dispersity,r.DH_pct,r.eta_mL_g,r.kH,r.kH_from_Kraemer,f.mass,f.dh,f.dh10Factor,f.trainingGrades.join('|'),r.source_doi].join(',');
    });
    const href=URL.createObjectURL(new Blob([[header,...lines].join('\n')+'\n'],{type:'text/csv;charset=utf-8'})),a=document.createElement('a');
    a.href=href;a.download='PVOH_hidroliz_numune_sinamasi.csv';a.click();setTimeout(()=>URL.revokeObjectURL(href),1000);
  };
  const full=result.full.dh;
  el('coefficients').textContent='Beş numunenin tümüyle uyumda a = '+fmt(full.a)+', β = '+fmt(full.beta)+'. Model, mol kütlesi ve çözücü sabitken hidroliz derecesinin +10 yüzde puan değişimine [η] için '+fmt(full.dh10Factor)+' çarpanı atar. Bu değer tek deney serisindeki koşullu ilişkidir.';
  el('grade').onchange=()=>{observation();validation();};el('scope').onchange=validation;el('concentration').oninput=observation;
  el('app').hidden=false;observation();validation();
})();
