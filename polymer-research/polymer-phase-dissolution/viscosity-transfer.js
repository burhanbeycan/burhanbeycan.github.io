(() => {
  'use strict';
  const el=id=>document.getElementById(id);
  const fmt=(n,d=4)=>n.toLocaleString('tr-TR',{maximumSignificantDigits:d});
  const q=new URLSearchParams(location.search), core=window.ViscosityTransfer;
  if (!window.VIS_INTRINSIC || !core) {el('err').hidden=false;return;}
  const pairs=core.makePairs(window.VIS_INTRINSIC.rows);
  for(const p of pairs) el('sample').add(new Option(p.id,p.id));
  el('sample').value=pairs.some(p=>p.id===q.get('sample'))?q.get('sample'):'IK1500-1';
  if(q.get('solvent')==='benzene')el('direction').value='benzene';
  let output=[];
  function draw(){
    const from=el('direction').value,split=el('split').value;
    output=core.benchmark(pairs,from,split);
    const row=output.find(p=>p.id===el('sample').value);
    const max=Math.max(row.observed,row.unscaled,row.adjusted),bars=el('bars');bars.replaceChildren();
    for(const [label,value] of [['Hedefte ölçülen',row.observed],['Yalnız Rɡ³ ile',row.unscaled],['Düzeltilmiş',row.adjusted]]){
      const line=document.createElement('div'),caption=document.createElement('span'),track=document.createElement('span'),fill=document.createElement('span'),number=document.createElement('strong');
      line.className='vt-row';track.className='vt-track';fill.className='vt-fill';
      fill.style.width=String(100*value/max)+'%';track.append(fill);caption.textContent=label;number.textContent=fmt(value);
      line.append(caption,track,number);bars.append(line);
    }
    el('reading').textContent='Rɡ '+fmt(row.sourceRg)+' nm’den '+fmt(row.targetRg)+' nm’ye değişiyor. Hacim oranı '+fmt(row.radiusRatio**3)+', ölçülen [η] oranı '+fmt(row.observed/row.sourceEta)+'. Düzeltme katsayısı '+fmt(row.factor)+' ('+row.train.join(', ')+').';
    const mape=k=>output.reduce((v,r)=>v+Math.abs(r[k]),0)/output.length;
    el('scores').textContent='Altı ayrılmış hedef ölçüm için ortalama mutlak yüzde sapma yalnız Rɡ³ ile %'+fmt(mape('baselineError'),3)+', düzeltmeyle %'+fmt(mape('adjustedError'),3)+'.';
    const tbody=el('results');tbody.replaceChildren();
    for(const r of output){
      const tr=tbody.insertRow();if(r.id===row.id)tr.className='vt-picked';
      [r.id,fmt(r.observed),fmt(r.unscaled),fmt(r.adjusted),(r.adjustedError>0?'+':'')+fmt(r.adjustedError,3)+' %'].forEach(v=>tr.insertCell().textContent=v);
    }
    const next=new URLSearchParams({source:'pair',system:'ps',sample:row.id,solvent:from,T:from==='benzene'?'298.15':'307.65'});
    const pairQuery=new URLSearchParams({sample:row.id});
    if(q.has('c')&&q.get('c').trim()&&Number.isFinite(Number(q.get('c')))&&Number(q.get('c'))>=0){next.set('c',q.get('c'));pairQuery.set('c',q.get('c'));}
    el('visBack').href='viscosity.html?'+next;el('backPair').href='chain-solvents.html?'+pairQuery;
  }
  el('save').onclick=()=>{
    const keys=['id','parent','from','to','sourceRg','targetRg','sourceEta','observed','unscaled','adjusted','baselineError','adjustedError','factor','train'];
    const header='sample,parent,source_solvent,target_solvent,Rg_source_nm,Rg_target_nm,eta_source_mL_g,eta_measured_target_mL_g,eta_volume_mL_g,eta_corrected_mL_g,volume_error_percent,corrected_error_percent,correction_factor,training_ids';
    const lines=output.map(r=>keys.map(k=>Array.isArray(r[k])?r[k].join('|'):r[k]).join(','));
    const url=URL.createObjectURL(new Blob([[header,...lines].join('\n')+'\n'],{type:'text/csv;charset=utf-8'}));
    const a=document.createElement('a');a.href=url;a.download='PS_Rg3_viskozite_'+el('split').value+'.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  for(const id of ['sample','direction','split'])el(id).onchange=draw;
  el('transferApp').hidden=false;draw();
})();
