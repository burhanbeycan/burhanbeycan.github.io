'use strict';
(() => {
  const el=id=>document.getElementById(id);
  const ctx=el('field').getContext('2d');
  const presets={quench:{N:100,chi:.85,phi:.15},slab:{N:100,chi:.30,phi:.20}};
  const params=new URLSearchParams(location.search);
  let linked=['N','chi','chi1','phi'].some(k=>params.has(k));
  const launchMode=params.get('mode')==='slab'?'slab':'quench';
  const linkedSpec={...presets[launchMode]};
  for(const k of ['N','chi','chi1','phi']) if(params.has(k)) linkedSpec[k]=Number(params.get(k));
  const temperature=params.has('T')?Number(params.get('T')):null;
  const thermal=temperature!==null && Number.isFinite(temperature) && temperature>0;
  const number=(v,max=3)=>v.toLocaleString('tr-TR',{maximumFractionDigits:max});
  const sci=v=>Math.abs(v)<1e-16?'≈ 0':v.toExponential(2).replace('.',',');
  const colors=[[20,62,83],[35,117,120],[148,194,148],[240,221,160],[176,101,52]];
  const palette=Array.from({length:256},(_,i)=>{
    const t=i/255*4,j=Math.min(3,Math.floor(t));
    return colors[j].map((v,k)=>Math.round(v+(colors[j+1][k]-v)*(t-j)));
  });
  let sim,mode=launchMode,active=false,request=0,tickAt=0,sourceTemperature=thermal,states={};
  function pause(){
    active=false;cancelAnimationFrame(request);
    el('play').textContent='▶ Oynat';el('play').setAttribute('aria-pressed','false');
  }
  function make(spec,nextMode,temp=false){
    const next=PolyCH.setup({...spec,mode:nextMode,dt:nextMode==='slab'?.12:.5});
    pause();sim=next;mode=nextMode;sourceTemperature=temp;
    states[mode]={spec:{N:sim.N,chi:sim.chi,chi1:sim.chi1,phi:sim.phi},thermal:temp};
    el('quench').setAttribute('aria-pressed',String(mode==='quench'));
    el('slab').setAttribute('aria-pressed',String(mode==='slab'));
    el('n').value=sim.N;el('x').value=sim.chi;el('c').value=sim.chi1;el('fraction').value=sim.phi;el('exact').value=sim.phi;
    el('warning').textContent='';render();
  }
  function graph(xs,ys,xLabel,yLabel,yMin,yMax,xMin,xMax,guide=null){
    const svg=el('plot');
    const tx=x=>50+(x-xMin)/(xMax-xMin)*415,ty=y=>138-(y-yMin)/(yMax-yMin)*119;
    const grid=[0,.5,1].map(k=>{
      const y=yMin+(yMax-yMin)*k;
      return '<path d="M50 '+ty(y)+'h415" stroke="#dbe4da"/><text x="43" y="'+(ty(y)+4)+'" text-anchor="end">'+number(y,2)+'</text>';
    }).join('');
    const xt=[0,.5,1].map(k=>{const x=xMin+(xMax-xMin)*k;return '<text x="'+tx(x)+'" y="157" text-anchor="middle">'+number(x,2)+'</text>';}).join('');
    const path=xs.map((x,i)=>(i?'L':'M')+tx(x).toFixed(2)+' '+ty(ys[i]).toFixed(2)).join(' ');
    const marker=guide!==null&&guide>=xMin&&guide<=xMax?'<path d="M'+tx(guide)+' 20 V138" stroke="#b06534" stroke-width="1.8" stroke-dasharray="5 4"/>':'';
    svg.innerHTML='<g font-family="system-ui" font-size="11" fill="#416054">'+grid+xt+'<text x="50" y="13">'+yLabel+'</text><text x="465" y="171" text-anchor="end">'+xLabel+'</text></g>'+marker+'<path d="'+path+'" fill="none" stroke="#1b7560" stroke-width="2.3"/>';
  }
  function rows(entries) {
    const dl=el('measure');dl.replaceChildren();
    for(const [name,value] of entries){
      const block=document.createElement('div'), key=document.createElement('dt'), val=document.createElement('dd');
      key.textContent=name;val.textContent=value;block.append(key,val);dl.append(block);
    }
  }
  function spectrum(){
    const n=sim.size,maxBin=n/2,bins=new Float64Array(maxBin+1),counts=new Int32Array(maxBin+1);
    for(let i=1;i<sim.real.length;i++){
      const k=Math.round(Math.sqrt(sim.shell[i]));
      if(k>maxBin) continue;
      bins[k]+=sim.real[i]**2+sim.imag[i]**2;counts[k]++;
    }
    const xs=[],ys=[],area=n*n;
    for(let j=1;j<bins.length;j++) if(counts[j]&&bins[j]>0) {xs.push(j*2*Math.PI/n);ys.push(Math.log10(bins[j]/(area*counts[j])));}
    const ymin=ys.length?Math.max(-14,Math.floor(Math.min(...ys))):-14;
    const ymax=ys.length?Math.max(ymin+1,Math.ceil(Math.max(...ys))):0;
    graph(xs,ys.map(y=>Math.max(ymin,y)),'q* (ters uzunluk)','log₁₀ S(q*)',ymin,ymax,0,Math.PI,sim.currentPrediction?.q??null);
    el('plot').setAttribute('aria-label','Dalga sayısına göre radyal ortalama yapı faktörü');
    el('graphTitle').textContent='Dalgalanmaların dalga boyu';
    el('graphNote').textContent=sim.currentPrediction?'Kesikli işaret, doğrusal teoride ızgaranın en hızlı büyüyen dalga sayısıdır. Çizgi, dalga sayısı bantlarında ortalanan güncel gücü gösterir.':'Çizgi dalga sayısı bantlarında ortalanan güncel gücü gösterir.';
    const spread=Math.sqrt(sim.field.reduce((sum,v)=>sum+(v-sim.phi)**2,0)/area);
    rows([['Bileşim dalgalanması · σφ',number(spread,7)],['Yerel eğrilik · f″',number(sim.curvature,5)],
      ['Beklenen dalga boyu · λ*',sim.currentPrediction?number(sim.currentPrediction.lambda,2):'Büyüyen mod yok'],
      ['Erken büyümeden λ*',sim.early?number(sim.early.lambda,2)+' (t* '+number(sim.early.t,2)+')':'—']]);
    el('interpret').textContent=sim.curvature>=0?'Küçük dalgalanmalar sönümlenir. Denge faz sayısını faz diyagramında kontrol edin.':sim.currentPrediction?'Seçili bileşim spinodal bölgede. Pozitif büyüme hızına sahip modlar güçlenirken kısa dalgalar sönümlenir.':'Yerel eğrilik negatif, ancak bu kutuda izinli dalga boyları arasında büyüyen mod bulunmuyor.';
  }
  function solvent(){
    const p=PolyCH.profile(sim), xs=Array.from(p,(_,x)=>x), ys=Array.from(p,v=>1-v);
    const metrics=PolyCH.slabMetrics(sim,p), initial=sim.initialSlab;
    graph(xs,ys,'x* (konum)','1 − φ',0,1,0,sim.size-1,metrics.x20);
    el('plot').setAttribute('aria-label','Tabaka derinliğine göre y yönünde ortalanmış çözücü hacim kesri');
    el('graphTitle').textContent='Tabaka boyunca çözücü profili';
    el('graphNote').textContent='Çizgi her konumdaki ortalama çözücü hacim kesridir. Kesikli işaret, çözücü kontrastının %20 eşiğindeki sol cepheyi izler.';
    rows([['Başlangıç tabakasında çözücü',number(metrics.solvent,5)],['Net çözücü alımı · Δ⟨1−φ⟩',number(metrics.solvent-initial.solvent,5)],
      ['Cephenin içe ilerlemesi · Δx*',metrics.x20===null?'Eşik seçilemiyor':number(metrics.x20-initial.x20,3)],
      ['Ara yüz genişliği · Δx*',metrics.width===null?'Eşik seçilemiyor':number(metrics.width,3)]]);
    el('interpret').textContent='Pozitif net alım, başlangıçta tabakanın bulunduğu bölgeye çözücü girdiğini gösterir. Tekil zincirin çözeltide serbest kalması ayrı bir gözlenebilirdir.';
  }
  function render(){
    if(!sim)return;
    let lo=Infinity, hi=-Infinity;
    for(const p of sim.field){lo=Math.min(lo,p);hi=Math.max(hi,p);}
    if(!el('contrast').checked){lo=0;hi=1;}
    if(hi-lo<1e-10){lo-=5e-11;hi+=5e-11;}
    const img=ctx.createImageData(sim.size,sim.size);
    for(let i=0;i<sim.field.length;i++){const rgb=palette[Math.max(0,Math.min(255,Math.round((sim.field[i]-lo)*255/(hi-lo))))];for(let k=0;k<3;k++)img.data[4*i+k]=rgb[k];img.data[4*i+3]=255;}
    ctx.putImageData(img,0,0);
    const digits=hi-lo<.001?7:5;
    el('lo').textContent=number(lo,digits);el('hi').textContent=number(hi,digits);el('scaleTitle').textContent='Yerel φ';
    el('time').textContent='t* = '+number(sim.time,2);
    el('fractionValue').textContent=number(sim.phi,6);
    el('cond').textContent='N = '+number(sim.N,2)+(sim.chi1?' · χ(φ) = '+number(sim.chi,5)+(sim.chi1<0?' − ':' + ')+number(Math.abs(sim.chi1),5)+'φ':' · χ = '+number(sim.chi,5))+(sourceTemperature?' · Aktarılan sıcaklık '+number(temperature,2)+' K':' · Boyutsuz model koşulu');
    el('mass').textContent=sci(PolyCH.mean(sim.field)-sim.initialMean);el('energy').textContent=sci(sim.energy-sim.initialEnergy);
    if(mode==='quench'){
      el('initialInfo').textContent='Yaklaşık homojen başlangıçta küçük, sabit tohumlu bileşim dalgalanmaları var. Renk uçları alandaki gerçek φ değerleridir.';spectrum();
    }else{
      el('initialInfo').textContent='Ortadaki polimerce zengin tabaka, iki yüzünden çözücüce zengin ortama açılır. Ölçülen zaman ve mesafe indirgenmiş birimlerdedir.';solvent();
    }
  }
  function take(steps){
    try {for(let i=0;i<steps;i++)PolyCH.step(sim);render();}
    catch(e){pause();el('warning').textContent=e.message;}
  }
  function animate(t){
    if(!active)return;
    if(t-tickAt>=48){tickAt=t;take(3);}
    if(active)request=requestAnimationFrame(animate);
  }
  el('play').onclick=()=>{
    if(active)return pause();
    active=true;el('play').textContent='Ⅱ Duraklat';el('play').setAttribute('aria-pressed','true');tickAt=0;request=requestAnimationFrame(animate);
  };
  for(const m of ['quench','slab'])el(m).onclick=()=>{
    if(m===mode)return;
    const next=states[m]??(linked?{spec:linkedSpec,thermal}:{spec:presets[m],thermal:false});
    try {make(next.spec,m,next.thermal);}catch(e){el('warning').textContent=e.message;}
  };
  el('step').onclick=()=>{pause();take(1);};
  el('reset').onclick=()=>make(states[mode].spec,mode,states[mode].thermal);
  el('contrast').onchange=render;
  el('fraction').oninput=()=>{try{make({...states[mode].spec,phi:+el('fraction').value},mode,states[mode].thermal);}catch(e){el('warning').textContent=e.message;}};
  el('apply').onclick=()=>{
    if(!['n','x','c','exact'].every(id=>el(id).reportValidity()))return;
    try{const n=+el('n').value,x=+el('x').value,c=+el('c').value;make({N:n,chi:x,chi1:c,phi:+el('exact').value},mode,states[mode].thermal && x===sim.chi && c===sim.chi1 && n===sim.N);}catch(e){el('warning').textContent=e.message;}
  };
  document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});
  try{make(linkedSpec,mode,thermal);}catch(e){linked=false;make(presets[mode],mode);el('warning').textContent='Aktarılan değer çalışma aralığı dışında. Ekranda boyutsuz başlangıç örneği açık.';}
})();
