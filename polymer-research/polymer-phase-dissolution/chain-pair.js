(() => {
  'use strict';
  const $=id=>document.getElementById(id), f=(n,k=4)=>n.toLocaleString('tr-TR',{maximumSignificantDigits:k});
  if(!window.ChainPhysics || !window.CHAIN_DATA || !window.PS_BENZENE_25C){$('pairError').hidden=false;return;}
  const pairs=window.PS_BENZENE_25C.map(benzene=>({benzene,theta:window.CHAIN_DATA.rows.find(r=>r.id===benzene.id)}));
  function stats(pair,c) {
    const under=r=>window.ChainPhysics.overlap(r.Mw_g_mol,r.Rg_nm,c);
    const alpha=pair.benzene.Rg_nm/pair.theta.Rg_nm;
    return {theta:under(pair.theta),benzene:under(pair.benzene),alpha};
  }
  pairs.forEach(p=>$('pairSample').add(new Option(p.theta.id+' · '+f(p.theta.Mw_g_mol/1e6)+' MDa (θ)',p.theta.id)));
  const params=new URLSearchParams(location.search);
  $('pairSample').value=pairs.some(p=>p.theta.id===params.get('sample'))?params.get('sample'):'IK1500-1';
  if(params.has('c') && params.get('c').trim()!=='' && Number.isFinite(+params.get('c')) && +params.get('c')>=0) $('pairC').value=String(+params.get('c'));
  function plotRadius(circle,line,cx,r) {
    const cy=197-r;
    for(const [a,n] of Object.entries({cx,cy,r}))circle.setAttribute(a,n);
    for(const [a,n] of Object.entries({x1:cx,y1:cy,x2:cx+r,y2:cy}))line.setAttribute(a,n);
  }
  function update() {
    const c=$('pairC').valueAsNumber, ok=Number.isFinite(c) && c>=0;
    $('pairInvalid').hidden=ok;$('pairResults').hidden=!ok;if(!ok)return;
    const pair=pairs.find(p=>p.theta.id===$('pairSample').value),s=stats(pair,c);
    const scale=90/506;
    plotRadius($('ringC'),$('radC'),130,pair.theta.Rg_nm*scale);
    plotRadius($('ringB'),$('radB'),365,pair.benzene.Rg_nm*scale);
    $('alpha').textContent=f(s.alpha)+' kat';$('volumeRatio').textContent=f(s.alpha**3)+' kat';
    const status=(r)=>r===1?'eşikle aynı':r>1?'eşiğin üzerinde':'eşiğin altında';
    $('reading').textContent=f(c)+' g/L’de, seçilen geometrik tanıma göre siklohekzandaki çözelti '+status(s.theta.ratio)+', benzendeki çözelti '+status(s.benzene.ratio)+'.';
    const body=$('pairBody');body.replaceChildren();
    for(const [key,label] of [['theta','Siklohekzan · 34,5 °C'],['benzene','Benzen · 25 °C']]){
      const row=body.insertRow();
      for(const v of [label,f(pair[key].Mw_g_mol/1e6,4),f(pair[key].Rg_nm),f(s[key].cstar_g_L,4),f(s[key].ratio,4)])row.insertCell().textContent=v;
    }
    const back='chain.html?'+new URLSearchParams({source:'pair',system:'ps',sample:pair.theta.id,T:'307.65',c:String(c)});
    $('chainTop').href=back;$('chainFromPair').href=back;
  }
  $('pairSample').onchange=update;$('pairC').oninput=update;
  $('pairCsv').onclick=()=>{
    const c=$('pairC').valueAsNumber;if(!Number.isFinite(c)||c<0)return;
    const lines=['sample,source_doi,solvent,T_K,Mw_g_mol,Rgz_nm,cstar_g_L,c_g_L,c_over_cstar,Rg_benzene_over_Rg_theta'];
    for(const p of pairs) {const s=stats(p,c);for(const key of ['theta','benzene'])lines.push([p.theta.id,'10.1021/ma60066a022',key==='theta'?'cyclohexane':'benzene',key==='theta'?307.65:298.15,p[key].Mw_g_mol,p[key].Rg_nm,s[key].cstar_g_L,c,s[key].ratio,s.alpha].join(','));}
    const url=URL.createObjectURL(new Blob([lines.join('\n')+'\n'],{type:'text/csv;charset=utf-8'}));
    const a=document.createElement('a');a.href=url;a.download='PS_iki_cozucu_alti_fraksiyon.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  $('pairApp').hidden=false;update();
})();
