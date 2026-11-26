(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const model = window.SolventActivity;
  if (!model || !window.DenseFH || !window.DENSE_RESULTS) {
    $('error').hidden = false;
    return;
  }
  const specs = {
    activity: ['Siklohekzanın aktivitesi', 'aₛ', ''],
    minusLogA: ['Aktivitenin logaritmik azalması', '−ln aₛ', ''],
    deltaMu_J_mol: ['Çözücü kimyasal potansiyelindeki değişim', 'Δμₛ (J/mol)', ' J/mol']
  };
  const f = (x, digits=5) => x.toLocaleString('tr-TR', {maximumSignificantDigits:digits});
  const rows = model.compare();
  const grids = new Map();
  const svgNS = 'http://www.w3.org/2000/svg';
  const svg = $('isotherm');
  const make = (tag, attrs, txt) => {
    const el = document.createElementNS(svgNS, tag);
    for (const [k,v] of Object.entries(attrs)) el.setAttribute(k, v);
    if (txt !== undefined) el.textContent = txt;
    svg.append(el); return el;
  };
  const points = T => {
    if (!grids.has(T)) {
      grids.set(T, Array.from({length:141}, (_,i) => {
        const phi = i/200;
        return {phi, baseline:model.predict(T,phi,'baseline'), extended:model.predict(T,phi,'extended')};
      }));
    }
    return grids.get(T);
  };
  function render() {
    const T = Number($('activityT').value), key = $('quantity').value, p = Number($('phi').value);
    const [title, axis, unit] = specs[key], series = rows.filter(r=>r.T_K===T);
    const sweep = points(T);
    const old = model.predict(T,p,'baseline')[key], ext = model.predict(T,p,'extended')[key];
    const match = series.find(r=>Math.abs(r.phi-p)<1e-9);
    const boundValues = [...sweep.flatMap(r=>[r.extended[key],r.baseline[key]]), ...series.map(r=>r.observed[key])];
    const min = Math.min(...boundValues), max=Math.max(...boundValues), padding=(max-min)*.035;
    const low=min-padding, high=max+padding;
    const left=73,right=575,top=22,bottom=248, X=x=>left+(right-left)*x/.7, Y=y=>bottom-(bottom-top)*(y-low)/(high-low);
    svg.replaceChildren();
    const text = (x,y,s,attrs={})=>make('text',{x,y,fill:'#466255','font-family':'system-ui,sans-serif','font-size':14,...attrs},s);
    // Bounds come from the pressure fit, g/L × mL/g / 1000.
    const loPhi=8.1*.934/1000, hiPhi=259*.934/1000;
    make('rect',{x:X(loPhi),y:top,width:X(hiPhi)-X(loPhi),height:bottom-top,fill:'#e6f1e9'});
    for(let j=0;j<=4;j++) {
      const y=min+j*(max-min)/4;
      make('line',{x1:left,y1:Y(y),x2:right,y2:Y(y),stroke:'#d2dcd4'});
      text(left-8,Y(y)+5,f(y,4),{'text-anchor':'end'});
    }
    for (let i=0;i<=7;i++) text(X(i/10),bottom+21,f(i/10),{'text-anchor':'middle'});
    text((left+right)/2,293,'Polimer hacim kesri · φ',{'text-anchor':'middle'});
    text(left,15,axis,{'font-weight':'600'});
    make('line',{x1:X(p),x2:X(p),y1:top,y2:bottom,stroke:'#8a6955','stroke-dasharray':'2 4'});
    for (const [id,stroke,dash] of [['baseline','#416580','7 4'],['extended','#176b58','none']]) {
      const d=sweep.map((r,i)=>(i?'L':'M')+X(r.phi).toFixed(2)+' '+Y(r[id][key]).toFixed(2)).join(' ');
      make('path',{d,fill:'none',stroke,'stroke-dasharray':dash,'stroke-width':2.6});
    }
    for (const r of series) {
      const dot=make('circle',{cx:X(r.phi),cy:Y(r.observed[key]),r:match===r?6:4.5,fill:'white',stroke:'#984726','stroke-width':2.4});
      const t=document.createElementNS(svgNS,'title'); t.textContent='φ = '+f(r.phi)+' · '+f(r.observed[key])+unit; dot.append(t);
    }
    $('chartTitle').textContent=title;
    $('phiText').textContent='φ = '+f(p);
    $('condition').textContent=f(T-273.15,4)+' °C · φ = '+f(p);
    $('oldValue').textContent=f(old)+unit;
    $('newValue').textContent=f(ext)+unit;
    $('observedValue').classList.toggle('no',!match);
    $('observedValue').textContent=match?f(match.observed[key])+unit:'Bu bileşimde ölçüm yok';
    $('osmoticLink').href='osmotic.html?'+new URLSearchParams({sample:'ps-kg-fr2',model:'joint',T:String(T),phi:String(p)});
    const body=$('recordBody'); body.replaceChildren();
    for (const r of series) {
      const tr=body.insertRow();
      if (r===match) tr.className='selected';
      const btn=document.createElement('button');btn.type='button';btn.textContent=f(r.phi);
      btn.setAttribute('aria-label','φ = '+f(r.phi)+' bileşimini seç');
      btn.onclick=()=>{ $('phi').value=String(r.phi); render(); $('phi').focus(); };
      tr.insertCell().append(btn);
      for(const name of ['observed','baseline','extended']) tr.insertCell().textContent=f(r[name][key])+unit;
    }
    const errs=k=>series.reduce((sum,r)=>sum+Math.abs(r[k][key]-r.observed[key]),0)/series.length;
    $('assessment').textContent='Bu sıcaklıktaki '+series.length+' damıtma noktasında '+axis+' için ortalama mutlak hata önceki modelde '+f(errs('baseline'),3)+unit+', basınca uydurulan modelde '+f(errs('extended'),3)+unit+'.';
  }
  const params=new URLSearchParams(location.search);
  const T=Number(params.get('T'));
  if ([307.15,317.15].includes(T)) $('activityT').value=String(T);
  const p=Number(params.get('phi'));
  if (params.has('phi') && Number.isFinite(p) && p>=0 && p<=.7) $('phi').value=String(p);
  $('activityT').onchange=()=>{const available=rows.filter(r=>r.T_K===Number($('activityT').value));$('phi').value=String(available[0].phi);render();};
  $('quantity').onchange=render;
  $('phi').oninput=render;
  $('saveCsv').onclick=()=>{
    const csv=[['source_doi','paper_table','sample','T_K','Mn_g_mol','phi','observation_reported_minus_ln_a','observed_a','baseline_a','pressure_fit_a','observed_deltaMu_J_mol','baseline_deltaMu_J_mol','pressure_fit_deltaMu_J_mol'].join(',')];
    rows.forEach(r=>csv.push(['10.1021/ja01517a022','III','II',r.T_K,72000,r.phi,r.minusLogA,r.observed.activity,r.baseline.activity,r.extended.activity,r.observed.deltaMu_J_mol,r.baseline.deltaMu_J_mol,r.extended.deltaMu_J_mol].join(',')));
    const url=URL.createObjectURL(new Blob([csv.join('\n')+'\n'],{type:'text/csv;charset=utf-8'}));
    const a=document.createElement('a');a.href=url;a.download='PS_siklohekzan_aktivite_Tablo_III.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  $('work').hidden=false;
  render();
})();
