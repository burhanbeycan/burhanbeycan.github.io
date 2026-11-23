'use strict';

const $ = (id) => document.getElementById(id);
const S = {samples: new Map(), template: null, apps: [], active: null, material: null, inputs: {}, results: {}, runs: [], dirty: {}, busy: false};
const symbols = {chain:'⌁', transport:'⇥', network:'◇', electrolyte:'±', chemistry:'⬡', synthesis:'⇌'};
const shortNames = {chain:'Chain physics',transport:'Molecular transport',network:'Networks & swelling',electrolyte:'Ion transport',chemistry:'Molecular chemistry',synthesis:'Polymer synthesis'};
const colors = ['#167f86','#ba7c42','#7561a3','#477aaf','#b35e6b','#758b45'];
const disabledBeforeRun = new WeakMap();
const clone = (x) => JSON.parse(JSON.stringify(x));
const el = (tag, cls, text) => {const n=document.createElement(tag); if(cls)n.className=cls; if(text!==undefined)n.textContent=text; return n;};

async function api(path, data) {
  const response=await fetch(path, data===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});
  let value;
  try {value=await response.json();} catch {throw new Error('The local calculation service returned an unreadable response.');}
  if(!response.ok)throw new Error(typeof value.detail==='string'?value.detail:JSON.stringify(value.detail||value));
  return value;
}
function notify(text,error=false){$('notice').hidden=!text;$('notice').textContent=text||'';$('notice').classList.toggle('error',error);}
function toast(text){$('toast').textContent=text;$('toast').hidden=false;setTimeout(()=>{$('toast').hidden=true;},3500);}
function fmt(v){if(v===null||v===undefined)return '—';if(typeof v!=='number')return String(v);if(!Number.isFinite(v))return '—';if(v===0)return '0';const a=Math.abs(v);if(a>=1e6||a<0.001)return v.toExponential(3);return new Intl.NumberFormat('en-US',{maximumSignificantDigits:5}).format(v);}
function slug(s){return String(s).replace(/[^a-zA-Z0-9._-]+/g,'-').slice(0,100);}
function download(name,content,type='application/json'){const blob=new Blob([content],{type});const a=el('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),2000);}
function exportJSON(name,doc){download(name,JSON.stringify(doc,null,2));}
function current(){return S.apps.find(a=>a.id===S.active);}
function markAllDirty(){for(const app of S.apps)S.dirty[app.id]=true;}
function lockUI(locked){for(const node of document.querySelectorAll('button,input,select,textarea')){if(locked){disabledBeforeRun.set(node,node.disabled);node.disabled=true;}else if(disabledBeforeRun.has(node)){node.disabled=disabledBeforeRun.get(node);disabledBeforeRun.delete(node);}}}

function setView(view){
  if(S.busy)return;
  for(const n of document.querySelectorAll('.view'))n.hidden=n.id!==`view-${view}`;
  for(const b of document.querySelectorAll('[data-view]'))b.setAttribute('aria-selected',String(b.dataset.view===view));
  if(view==='science')renderScience();
  if(view==='pilots')renderPilots();
  if(view==='material')renderMaterial();
  if(view==='workflow')renderWorkflow();
}

function updateSample(){
  S.samples.set(S.material.sample_id, clone(S.material));
  $('sample-name').textContent=S.material.name;
  $('sample-id').textContent=S.material.sample_id;
  $('run-count').textContent=String(S.runs.length);
}

function selectApp(id, autoRun=false){
  if(S.busy)return;
  const app=S.apps.find(a=>a.id===id);if(!app)return;
  S.active=id;
  document.documentElement.style.setProperty('--accent',app.accent||'#177b83');
  document.title=`${app.name} · Research Studio`;
  $('app-title').textContent=app.name;
  $('app-description').textContent=app.description;
  $('eyebrow').textContent=app.tagline.toUpperCase();
  $('breadcrumb').textContent=shortNames[id].toUpperCase();
  for(const b of $('app-nav').children)b.classList.toggle('active',b.dataset.app===id);
  if(!S.inputs[id])S.inputs[id]=Object.fromEntries(app.fields.map(f=>[f.key,f.default]));
  $('preset').replaceChildren(el('option',null,'Custom inputs'));
  $('preset').firstChild.value='';
  app.presets.forEach((p,i)=>{const n=el('option',null,p.name);n.value=String(i);$('preset').append(n);});
  $('preset-description').textContent='Initial values are demonstrations. Review coefficients for your own material.';
  renderFields();
  setView('analysis');
  notify('');
  if(S.results[id])renderResults(S.results[id]);else renderEmpty();
  $('run-state').textContent=S.results[id]?(S.dirty[id]?'Inputs or material changed · recalculate':'Stored calculation · inputs unchanged'):'Ready';
  updateExportState();
  if(autoRun&&!S.results[id])run();
}

function renderFields(){
  const app=current(), inputs=S.inputs[S.active], groups=new Map();
  $('fields').replaceChildren();
  const byKey=new Map(app.fields.map(f=>[f.key,f]));
  function visible(field,visited=new Set()){
    if(!field.show_if)return true;
    if(visited.has(field.key))return false;
    visited.add(field.key);
    const {key,value}=field.show_if, parent=byKey.get(key);
    const matches=Array.isArray(value)?value.includes(inputs[key]):value===inputs[key];
    return matches && Boolean(parent) && visible(parent,visited);
  }
  for(const f of app.fields){
    if(!visible(f))continue;
    const group=f.group||'Parameters';
    if(!groups.has(group)){const container=el('section','field-group');container.append(el('h3',null,group));groups.set(group,container);$('fields').append(container);}
    const row=el('div','field');
    const label=el('label',null,f.label);label.htmlFor=`field-${f.key}`;
    if(f.unit)label.append(el('span','unit',f.unit));
    let input;
    if(f.type==='select'){
      input=el('select');for(const option of f.options||[]){const n=el('option',null,option.label);n.value=option.value;input.append(n);}
    }else if(f.type==='textarea'){input=el('textarea');input.rows=5;}
    else {input=el('input');input.type=f.type==='text'?'text':'number';if(input.type==='number'){input.step=f.step??'any';if(f.min!==undefined)input.min=f.min;if(f.max!==undefined)input.max=f.max;}input.required=f.required??true;}
    input.id=`field-${f.key}`;input.name=f.key;input.value=inputs[f.key]??'';
    input.addEventListener(f.type==='select'?'change':'input',()=>{
      inputs[f.key]=f.type==='number'?(input.value===''?null:Number(input.value)):input.value;
      S.dirty[S.active]=true;updateExportState();$('run-state').textContent='Inputs changed · recalculate';
      if(f.type==='select')renderFields();
    });
    row.append(label,input);if(f.help){const help=el('p','field-help',f.help);help.id=`help-${f.key}`;input.setAttribute('aria-describedby',help.id);row.append(help);}
    groups.get(group).append(row);
  }
}

function renderEmpty(){
  const n=el('div','empty-state');n.append(el('div','empty-icon','∫'),el('h3',null,'A documented calculation starts here'),el('p',null,'Review the selected model and its assumptions, then run a calculation. The material record travels with each result.'));
  $('results').replaceChildren(n);
}
function updateExportState(){const disabled=!S.results[S.active]||Boolean(S.dirty[S.active])||S.busy;for(const id of ['export-result','export-csv','print-btn'])$(id).disabled=disabled;}

async function run(){
  if(S.busy)return;
  if(!$('model-form').reportValidity())return;
  const id=S.active;S.busy=true;lockUI(true);$('run-btn').firstChild.textContent='Computing… ';$('run-state').textContent='Solving the selected model';notify('');updateExportState();
  const start=performance.now();
  try{
    const doc=await api(`/api/calculate/${id}`,{inputs:clone(S.inputs[id]),material:S.material});
    S.results[id]=doc;S.runs.push(doc);if(S.runs.length>50)S.runs.shift();S.material=doc.material;markAllDirty();S.dirty[id]=false;
    updateSample();
    if(S.active===id){renderResults(doc);$('run-state').textContent=`Completed · ${((performance.now()-start)/1000).toFixed(2)} s`;}
  }catch(e){notify(e.message,true);$('run-state').textContent='Review model inputs';}
  finally{S.busy=false;lockUI(false);$('run-btn').disabled=false;$('run-btn').firstChild.textContent='Run calculation ';updateExportState();}
}

function renderResults(doc){
  const out=doc.results, root=$('results');root.replaceChildren();
  const metrics=el('div','metrics');
  for(const m of out.metrics||[]){const card=el('div','metric');const value=el('div','metric-value',fmt(m.value));if(m.unit)value.append(el('span','unit',m.unit));card.append(el('div','metric-label',m.label),value);if(m.note)card.append(el('div','metric-note',m.note));metrics.append(card);}
  root.append(metrics);
  if((out.warnings||[]).length){const card=el('div','card warning-card');card.append(el('h3',null,'Interpretation & scope'));const list=el('ul');out.warnings.forEach(w=>list.append(el('li',null,w)));card.append(list);root.append(card);}
  const grouped=new Map();
  for(const view of out.molecular_views||[]){
    if(view.comparison_group){
      if(!grouped.has(view.comparison_group)){const grid=el('div','molecular-comparison');grouped.set(view.comparison_group,grid);root.append(grid);}
      grouped.get(view.comparison_group).append(moleculeCard(view));
    }else root.append(moleculeCard(view));
  }
  for(const plot of out.plots||[])root.append(plotCard(plot));
  for(const table of out.tables||[])root.append(tableCard(table));
  const notes=el('div','card assumptions');notes.append(el('span','model-badge','MODEL CALCULATION'),el('h3',null,'Assumptions & scientific interpretation'));
  const list=el('ul');for(const note of out.notes||[])list.append(el('li',null,note));notes.append(list);
  const handoffs=el('div','handoff-actions');
  for(const h of out.handoffs||[]){const available=S.apps.some(a=>a.id===h.target);const button=el('button','primary small',h.label||`Transfer to ${h.target}`);button.disabled=!available;button.title=available?'Review the mapped inputs in the next application':'Install the target application in this workbench, or export and import the result JSON.';button.addEventListener('click',()=>transfer(doc,h.target));handoffs.append(button);}
  if(handoffs.children.length)notes.append(handoffs);
  root.append(notes);
  if(doc.study?.declaration){
    const details=el('details','card');details.append(el('summary',null,'Declared model and experimental counterpart'));
    details.append(tableCard({title:'Specimen '+doc.study.specimen,columns:['Field','Declared value'],rows:Object.entries(doc.study.declaration)}));
    for(const text of doc.study.audit)details.append(el('p',null,text));
    root.append(details);
  }
  const report=el('details','card print-only');report.open=true;report.append(el('summary',null,'Calculation inputs and sources'));
  const pre=el('pre','equation',JSON.stringify(doc.inputs,null,2));report.append(pre);for(const ref of doc.references||[])report.append(el('p',null,`${ref.title} — ${ref.url}`));root.append(report);
  root.append(el('div','run-provenance',`Sample ${doc.material.sample_id} · ${doc.application.name} ${doc.application.version} · ${new Date(doc.created_utc).toLocaleString()}\nRun ${doc.run_id} · SHA-256 ${doc.checksum_sha256}`));
}

function svgEl(tag, attrs={}, text){const n=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v]of Object.entries(attrs))n.setAttribute(k,String(v));if(text!==undefined)n.textContent=text;return n;}
function plotCard(plot){
  const card=el('div','chart-card'), heading=el('div','chart-heading');heading.append(el('h3',null,plot.title));const save=el('button','quiet','↓ SVG');heading.append(save);card.append(heading);
  const W=720,H=340,margin={l:80,r:22,t:25,b:65},pw=W-margin.l-margin.r,ph=H-margin.t-margin.b;
  const xlog=plot.x_scale==='log',ylog=plot.y_scale==='log';
  const valid=(x,y)=>Number.isFinite(x)&&Number.isFinite(y)&&(!xlog||x>0)&&(!ylog||y>0);
  const series=(plot.series||[]).map(s=>({name:s.name,points:(s.x||[]).map((x,i)=>[x,(s.y||[])[i]]).filter(([x,y])=>valid(x,y))}));
  const points=series.flatMap(s=>s.points),svg=svgEl('svg',{viewBox:`0 0 ${W} ${H}`,role:'img','aria-label':`${plot.title}. ${plot.x_label}; ${plot.y_label}.`});
  svg.append(svgEl('title',{},plot.title));svg.append(svgEl('rect',{x:0,y:0,width:W,height:H,fill:'#ffffff'}));
  if(!points.length){card.append(el('p','subtle','No finite points are available for this plot.'));save.disabled=true;return card;}
  const tx=x=>xlog?Math.log10(x):x,ty=y=>ylog?Math.log10(y):y;
  let xmin=Math.min(...points.map(p=>tx(p[0]))),xmax=Math.max(...points.map(p=>tx(p[0]))),ymin=Math.min(...points.map(p=>ty(p[1]))),ymax=Math.max(...points.map(p=>ty(p[1])));
  if(xmax===xmin){xmin-=.5;xmax+=.5;}if(ymax===ymin){const d=Math.abs(ymax)*.05||.5;ymin-=d;ymax+=d;}else{const d=(ymax-ymin)*.05;ymin-=d;ymax+=d;}
  const X=x=>margin.l+(tx(x)-xmin)/(xmax-xmin)*pw,Y=y=>margin.t+ph-(ty(y)-ymin)/(ymax-ymin)*ph;
  const tick=(v,log)=>log?fmt(10**v):fmt(v);
  for(let i=0;i<=5;i++){
    const x=margin.l+i*pw/5,y=margin.t+i*ph/5;
    svg.append(svgEl('line',{x1:margin.l,y1:y,x2:W-margin.r,y2:y,stroke:'#e7edf1','stroke-width':1}));
    svg.append(svgEl('text',{x:margin.l-12,y:y+4,'text-anchor':'end',fill:'#71828c','font-family':'Arial, sans-serif','font-size':10},tick(ymax-i*(ymax-ymin)/5,ylog)));
    svg.append(svgEl('line',{x1:x,y1:margin.t+ph,x2:x,y2:margin.t+ph+5,stroke:'#a4b4bd','stroke-width':1}));
    svg.append(svgEl('text',{x,y:margin.t+ph+22,'text-anchor':'middle',fill:'#71828c','font-family':'Arial, sans-serif','font-size':10},tick(xmin+i*(xmax-xmin)/5,xlog)));
  }
  svg.append(svgEl('line',{x1:margin.l,y1:margin.t,x2:margin.l,y2:margin.t+ph,stroke:'#9bafb9','stroke-width':1}));
  svg.append(svgEl('line',{x1:margin.l,y1:margin.t+ph,x2:W-margin.r,y2:margin.t+ph,stroke:'#9bafb9','stroke-width':1}));
  const legend=el('div','chart-legend');
  series.forEach((s,i)=>{
    const color=colors[i%colors.length];
    const path=svgEl('path',{d:s.points.map((p,j)=>`${j?'L':'M'}${X(p[0]).toFixed(3)},${Y(p[1]).toFixed(3)}`).join(' '),fill:'none',stroke:color,'stroke-width':2,'stroke-linejoin':'round'});path.append(svgEl('title',{},s.name));svg.append(path);
    if(s.points.length<35)for(const p of s.points)svg.append(svgEl('circle',{cx:X(p[0]),cy:Y(p[1]),r:2.7,fill:color}));
    const label=el('span'),swatch=el('i');swatch.style.background=color;label.append(swatch,document.createTextNode(s.name||`Series ${i+1}`));legend.append(label);
  });
  svg.append(svgEl('text',{x:margin.l+pw/2,y:H-14,'text-anchor':'middle',fill:'#405867','font-family':'Arial, sans-serif','font-size':11},plot.x_label||''));
  svg.append(svgEl('text',{transform:`translate(18 ${margin.t+ph/2}) rotate(-90)`,'text-anchor':'middle',fill:'#405867','font-family':'Arial, sans-serif','font-size':11},plot.y_label||''));
  const cursor=svgEl('g',{visibility:'hidden'}),line=svgEl('line',{y1:margin.t,y2:margin.t+ph,stroke:'#a3b6be','stroke-dasharray':'3 3'}),dot=svgEl('circle',{r:4,fill:colors[0],stroke:'white','stroke-width':2});cursor.append(line,dot);svg.append(cursor);
  const caption=el('p','chart-caption',`${xlog?'Logarithmic':'Linear'} x-axis · ${ylog?'logarithmic':'linear'} y-axis · hover to inspect values`);
  svg.addEventListener('mousemove',e=>{const rect=svg.getBoundingClientRect(),px=(e.clientX-rect.left)*W/rect.width;const s=series.find(s=>s.points.length);if(!s)return;let nearest=s.points[0];for(const p of s.points)if(Math.abs(X(p[0])-px)<Math.abs(X(nearest[0])-px))nearest=p;line.setAttribute('x1',X(nearest[0]));line.setAttribute('x2',X(nearest[0]));dot.setAttribute('cx',X(nearest[0]));dot.setAttribute('cy',Y(nearest[1]));cursor.setAttribute('visibility','visible');caption.textContent=`${s.name}: x = ${fmt(nearest[0])}, y = ${fmt(nearest[1])}`;});
  svg.addEventListener('mouseleave',()=>cursor.setAttribute('visibility','hidden'));
  save.addEventListener('click',()=>{const copy=svg.cloneNode(true);copy.lastChild.remove();download(`${slug(plot.title)}.svg`,new XMLSerializer().serializeToString(copy),'image/svg+xml');});
  card.append(svg,legend,caption);return card;
}

function moleculeCard(view){
  const card=el('div','chart-card molecule-card'), heading=el('div','chart-heading');
  heading.append(el('h3',null,view.title));const save=el('button','quiet','↓ Structure SVG');heading.append(save);card.append(heading);
  const caption=el('p','molecule-caption',view.caption||'Finite molecular connectivity diagram. Drawing coordinates are not an equilibrium conformation.');card.append(caption);
  const controls=el('div','molecule-controls');
  function toggle(label,initial){const wrap=el('label'),box=el('input');box.type='checkbox';box.checked=initial;wrap.append(box,document.createTextNode(label));controls.append(wrap);return box;}
  const allAtoms=toggle('Label carbon atoms',false), indices=toggle('Atom indices',false), units=toggle('Repeat-unit colours',true);
  const zoomLabel=el('label',null,'Zoom '),zoom=el('input');zoom.type='range';zoom.min='100';zoom.max='400';zoom.step='25';zoom.value='100';zoom.setAttribute('aria-label','Molecular drawing zoom');zoomLabel.append(zoom);controls.append(zoomLabel);card.append(controls);
  const viewport=el('div','molecule-viewport');viewport.tabIndex=0;viewport.setAttribute('aria-label','Molecular connectivity diagram. Scroll horizontally after zooming.');card.append(viewport);
  let svg;
  function draw(){
    const pad=45,atoms=view.atoms||[],bonds=view.bonds||[];
    const byIndex=new Map(atoms.map(a=>[a.index,a])), endColors={alpha:'#13766d',omega:'#a74679'};
    const xs=atoms.map(a=>a.x),ys=atoms.map(a=>a.y),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
    const W=Math.max(800,Math.min(12000,(xmax-xmin)*25+2*pad)),H=Math.max(360,Math.min(1600,(ymax-ymin)*25+2*pad));
    const scale=Math.min((W-2*pad)/Math.max(xmax-xmin,1),(H-2*pad)/Math.max(ymax-ymin,1));
    const X=a=>W/2+(a.x-(xmin+xmax)/2)*scale,Y=a=>H/2-(a.y-(ymin+ymax)/2)*scale;
    const font=Math.max(8,Math.min(17,scale*.52));
    svg=svgEl('svg',{viewBox:`0 0 ${W} ${H}`,role:'img','aria-label':`${view.title}. ${atoms.length} atoms, ${bonds.length} bonds.`});
    svg.style.width=view.comparison_group?`${zoom.value}%`:`${W*Number(zoom.value)/100}px`;svg.style.maxWidth='none';svg.append(svgEl('title',{},view.title),svgEl('desc',{},view.caption||'Finite molecular connectivity drawing.'),svgEl('rect',{width:W,height:H,fill:'#fff'}));
    if(units.checked)for(const atom of atoms){if(Number.isInteger(atom.unit_index)){svg.append(svgEl('circle',{cx:X(atom),cy:Y(atom),r:Math.max(5,font*.8),fill:colors[((atom.unit_index%colors.length)+colors.length)%colors.length],opacity:.11}));}}
    for(const atom of atoms){if(endColors[atom.end_role])svg.append(svgEl('circle',{cx:X(atom),cy:Y(atom),r:font*.95,fill:endColors[atom.end_role],opacity:.28}));}
    const line=(x1,y1,x2,y2,extra={})=>svg.append(svgEl('line',{x1,y1,x2,y2,stroke:'#2b4250','stroke-width':Math.max(1,Math.min(1.8,scale*.07)),'stroke-linecap':'round',...extra}));
    for(const bond of bonds){
      const a=byIndex.get(bond.begin),b=byIndex.get(bond.end);if(!a||!b)continue;
      const x1=X(a),y1=Y(a),x2=X(b),y2=Y(b),dx=x2-x1,dy=y2-y1,len=Math.hypot(dx,dy)||1,px=-dy/len,py=dx/len,sep=Math.max(1.1,Math.min(3,scale*.1));
      if(endColors[bond.end_role])line(x1,y1,x2,y2,{stroke:endColors[bond.end_role],'stroke-width':Math.max(4,font*.42),opacity:.55});
      if(bond.direction==='BEGINWEDGE')svg.append(svgEl('polygon',{points:`${x1},${y1} ${x2+px*sep*1.5},${y2+py*sep*1.5} ${x2-px*sep*1.5},${y2-py*sep*1.5}`,fill:'#2b4250'}));
      else if(bond.direction==='BEGINDASH'){for(let j=1;j<=7;j++){const t=j/8,cx=x1+t*dx,cy=y1+t*dy;line(cx-px*sep*t*1.5,cy-py*sep*t*1.5,cx+px*sep*t*1.5,cy+py*sep*t*1.5);}}
      else if(bond.order===2){for(const offset of [-sep/2,sep/2])line(x1+px*offset,y1+py*offset,x2+px*offset,y2+py*offset);}
      else if(bond.order===3){for(const offset of [-sep,0,sep])line(x1+px*offset,y1+py*offset,x2+px*offset,y2+py*offset);}
      else if(bond.order===1.5){line(x1,y1,x2,y2);line(x1+px*sep,y1+py*sep,x2+px*sep,y2+py*sep,{'stroke-dasharray':'3 3'});}
      else line(x1,y1,x2,y2);
    }
    const elementColors={O:'#be4d46',N:'#3b67af',S:'#9d831e',F:'#43916b',Cl:'#43916b',Br:'#8b5650',P:'#9a651b'};
    for(const atom of atoms){
      const group=svgEl('g');group.append(svgEl('title',{},`Atom ${atom.index}: ${atom.element}${Number.isInteger(atom.unit_index)?` · repeat ${atom.unit_index}`:''}${atom.chiral?` · ${atom.chiral}`:''}`));
      const chemicalLabel=atom.label||atom.element,show=atom.end_role||atom.element!=='C'||allAtoms.checked||atom.charge||atom.isotope;
      if(show){const text=svgEl('text',{x:X(atom),y:Y(atom),'text-anchor':'middle','dominant-baseline':'central','font-family':'Arial, sans-serif','font-size':font,'font-weight':500,fill:elementColors[atom.element]||'#2b4250',stroke:'#fff','stroke-width':font*.29,'paint-order':'stroke'},chemicalLabel);group.append(text);}
      else group.append(svgEl('circle',{cx:X(atom),cy:Y(atom),r:Math.max(3,font*.4),fill:'transparent'}));
      if(indices.checked)group.append(svgEl('text',{x:X(atom)+font*.4,y:Y(atom)-font*.65,'font-family':'Arial, sans-serif','font-size':font*.6,fill:'#71818b'},String(atom.index)));
      svg.append(group);
    }
    viewport.replaceChildren(svg);
  }
  for(const input of [allAtoms,indices,units,zoom])input.addEventListener('input',draw);
  draw();
  const legend=el('div','molecule-legend');legend.append(el('span','model-badge','FINITE MOLECULAR GRAPH'),el('span',null,'Lines are covalent bonds · colours identify repeat membership · hover to inspect atoms · scroll for large structures'));if(view.atoms.some(a=>a.end_role))legend.append(el('b','alpha-end','α initiation end'),el('b','omega-end','ω termination end'));card.append(legend);
  if(view.smiles){const details=el('details','molecule-smiles');details.append(el('summary',null,'Canonical isomeric SMILES'),el('p','equation',view.smiles));card.append(details);}
  save.addEventListener('click',()=>{const copy=svg.cloneNode(true);copy.removeAttribute('style');download(`${slug(view.title)}.svg`,new XMLSerializer().serializeToString(copy),'image/svg+xml');});
  return card;
}

function tableCard(table){
  const card=el('div','card');card.append(el('h3',null,table.title));const wrap=el('div','table-wrap'),t=el('table'),thead=el('thead'),tr=el('tr');
  for(const col of table.columns||[])tr.append(el('th',null,typeof col==='object'?col.label:col));thead.append(tr);t.append(thead);const body=el('tbody');
  for(const row of (table.rows||[]).slice(0,250)){const r=el('tr');const values=Array.isArray(row)?row:(table.columns||[]).map(c=>row[typeof c==='object'?c.key:c]);for(const value of values)r.append(el('td',null,typeof value==='object'&&value!==null?JSON.stringify(value):fmt(value)));body.append(r);}
  t.append(body);wrap.append(t);card.append(wrap);if((table.rows||[]).length>250)card.append(el('p','subtle','Showing 250 rows. The complete calculation record is available in the JSON export.'));return card;
}

function renderScience(){
  const app=current();if(!app)return;const root=$('view-science');root.replaceChildren();const grid=el('div','science-grid'),main=el('div'),side=el('div');
  const scope=el('div','card assumptions');scope.append(el('span','model-badge','SUPPORTED SCIENTIFIC SCOPE'),el('h3',null,app.name),el('p',null,app.scope));main.append(scope);
  for(const model of app.models||[]){const card=el('div','card');card.append(el('h3',null,model.name),el('div','equation',model.equation));const ul=el('ul');for(const a of model.assumptions||[])ul.append(el('li',null,a));card.append(ul);if(model.references?.length)card.append(el('p','subtle','Sources: '+model.references.join(' · ')));main.append(card);}
  const refs=el('div','card');refs.append(el('h3',null,'Primary sources & documentation'));
  for(const r of app.references||[]){const item=el('div','reference-item');item.append(el('span',null,r.id));const a=el('a',null,r.title);try{const url=new URL(r.url);if(['https:','http:'].includes(url.protocol))a.href=url.href;}catch{}a.target='_blank';a.rel='noopener noreferrer';item.append(a);refs.append(item);}side.append(refs);
  const verification=el('div','card');verification.append(el('h3',null,'How to interpret validation'),el('p',null,'The repository includes analytical, numerical and invalid-input checks. These establish behavior of the implemented equations. A named polymer requires independent parameter calibration and experimental validation.'),el('p',null,'Molecular Chemistry builds finite, valence-checked structures from explicit templates. Physical coefficients still require independent evidence; a structure drawing does not establish an equilibrium conformation or bulk property.'));side.append(verification);grid.append(main,side);root.append(grid);
}

function renderMaterial(){
  const root=$('view-material');root.replaceChildren();const grid=el('div','material-grid'),card=el('div','card');card.append(el('h3',null,'A polymer sample is a population'),el('p',null,'Keep chemistry, architecture, sample history and measurement provenance attached to every calculation. Unspecified information stays unspecified.'));
  const form=el('form','material-form');const fields=[['sample_id','Sample identifier'],['name','Sample name'],['chemistry.repeat_unit','Repeat unit / chemistry'],['chemistry.architecture','Architecture'],['chemistry.tacticity','Tacticity'],['chemistry.sequence','Sequence / composition'],['chemistry.end_groups','End groups'],['chemistry.repeat_smiles','Mapped repeat-unit SMILES'],['chemistry.stereochemistry_source','Stereochemistry evidence'],['chemistry.branch_topology','Branch / junction topology'],['chemistry.functional_groups','Functional groups'],['chemistry.synthesis_route','Polymerization route'],['chemistry.conversion_basis','Conversion / composition basis'],['conditions.temperature_K','Reference temperature (K)'],['conditions.state','Physical state'],['conditions.solvent','Solvent / environment'],['chemistry.notes','Chemical characterization notes'],['conditions.notes','Preparation & condition notes']];
  for(const [key,labelText] of fields){const parts=key.split('.'),value=parts.length===2?S.material[parts[0]][parts[1]]:S.material[key],row=el('div','field'+(key.endsWith('notes')?' wide':'')),label=el('label',null,labelText);label.htmlFor=`material-${key}`;const input=el(key.endsWith('notes')?'textarea':'input');input.id=`material-${key}`;input.name=key;if(key==='sample_id')input.readOnly=true;input.value=value??'';if(key.endsWith('temperature_K')){input.type='number';input.min='0.001';input.step='any';}else input.maxLength=key.endsWith('notes')?4000:200;row.append(label,input);form.append(row);}
  const save=el('button','primary small','Save material record');save.type='submit';form.append(save);form.addEventListener('submit',async e=>{e.preventDefault();const m=clone(S.material);for(const [key]of fields){const v=form.elements.namedItem(key).value,parts=key.split('.');const value=key.endsWith('temperature_K')?(v===''?null:Number(v)):v;if(parts.length===2)m[parts[0]][parts[1]]=value;else m[key]=value;}try{S.material=await api('/api/exchange/validate',m);updateSample();for(const a of S.apps)S.dirty[a.id]=true;updateExportState();toast('Material record saved. Recalculate to associate new results.');}catch(e){notify(e.message,true);}});card.append(form);
  const side=el('div'),scope=el('div','card assumptions');scope.append(el('h3',null,'What transfers between applications'),el('p',null,'The sample identity, chemistry metadata, population and calculation lineage travel together. Model coefficients retain their physical meaning and units. A tracer diffusion coefficient is not silently substituted for a Fickian transport coefficient.'));
  const exportButton=el('button','quiet','↓ Export material JSON');exportButton.addEventListener('click',()=>exportJSON(`${slug(S.material.sample_id)}-material.json`,S.material));scope.append(exportButton);side.append(scope);
  const dist=S.material.molar_mass_distribution;if(dist){side.append(tableCard({title:'Imported polymer population',columns:['M (g/mol)','Number fraction'],rows:dist.molar_mass_g_mol.map((m,i)=>[m,dist.number_fraction[i]])}));}
  const history=el('div','card');history.append(el('h3',null,'Calculation lineage'));if(!S.material.history.length)history.append(el('p',null,'No completed calculation is attached yet.'));else for(const h of S.material.history.slice(-10).reverse())history.append(el('p',null,`${h.application} · v${h.version} · ${h.run_id.slice(0,8)}`));side.append(history);grid.append(card,side);root.append(specimenPanel(),studyPanel(),evidencePanel(),grid);
}

function renderWorkflow(){
  const root=$('view-workflow');root.replaceChildren();root.append(el('p','workflow-intro','Each laboratory works independently. In this workbench, all installed laboratories share a versioned material record. Explicit handoffs map compatible quantities; general transfers carry material identity and provenance only.'));
  const strip=el('div','workflow-strip');
  for(const app of S.apps){const node=el('div','workflow-node');node.style.setProperty('--node-accent',app.accent);node.append(el('h3',null,shortNames[app.id]),el('p',null,app.tagline));const b=el('button','quiet','Open laboratory →');b.addEventListener('click',()=>selectApp(app.id,true));node.append(b);strip.append(node);}root.append(strip);
  const note=el('div','card assumptions');note.append(el('h3',null,'Defined transfers preserve chemical meaning'),el('p',null,'Chemistry Lab → Chain Lab: transfer the computed homopolymer repeat molar mass and explicit cap correction. The sample MWD, contour increment and stiffness require separate evidence. A finite oligomer is not a measured polymer population.'),el('p',null,'Network Lab → Transport Lab: transfer the fully swollen equilibrium film thickness and the user-supplied effective solute diffusivity. The receiving model treats a fixed geometry. Dynamic swelling is not implied by this handoff.'));root.append(note);
  root.append(el('h2',null,'Calculation history'));
  if(!S.runs.length){root.append(el('p','muted','Run an analysis to start a connected workflow.'));return;}
  for(const doc of [...S.runs].reverse()){
    const row=el('div','run-card'),info=el('div');info.append(el('h3',null,doc.application.name),el('p',null,`${doc.material.sample_id} · ${new Date(doc.created_utc).toLocaleString()}`),el('span','run-id',doc.run_id));const actions=el('div','run-actions');
    const open=el('button','quiet','View result');open.disabled=!S.apps.some(a=>a.id===doc.application.id);open.addEventListener('click',()=>{S.material=clone(S.samples.get(doc.material.sample_id)||doc.material);clearSpecimenInputs();markAllDirty();S.results[doc.application.id]=doc;S.inputs[doc.application.id]=clone(doc.inputs);S.dirty[doc.application.id]=false;updateSample();selectApp(doc.application.id);});actions.append(open);
    for(const target of S.apps.filter(a=>a.id!==doc.application.id)){const b=el('button','quiet',`→ ${shortNames[target.id]}`);b.addEventListener('click',()=>transfer(doc,target.id));actions.append(b);}
    const save=el('button','quiet','↓ JSON');save.addEventListener('click',()=>exportJSON(`${doc.application.id}-${doc.run_id.slice(0,8)}.json`,doc));actions.append(save);row.append(info,actions);root.append(row);
  }
}

async function transfer(doc,target){
  try{const result=await api('/api/exchange/transfer',{document:doc,target,destination:S.material});S.material=result.material;markAllDirty();selectApp(target);Object.assign(S.inputs[target],result.inputs);S.dirty[target]=true;renderFields();updateSample();updateExportState();notify(`Transferred from ${doc.application.name}. ${(result.notes||[]).join(' ')} Review the inputs, then run the receiving model.`);$('run-state').textContent='Transferred inputs · awaiting calculation';}
  catch(e){notify(e.message,true);}
}

function csvFor(doc){const rows=[['application',doc.application.name],['run_id',doc.run_id],['sample_id',doc.material.sample_id],[]];for(const view of doc.results.molecular_views||[]){rows.push(['molecular_view',view.title],['coordinates','arbitrary 2D drawing units, not molecular geometry'],['atom_index','element','drawing_x','drawing_y','unit_index','label','end_role']);for(const atom of view.atoms)rows.push([atom.index,atom.element,atom.x,atom.y,atom.unit_index??'',atom.label??'',atom.end_role??'']);rows.push([],['bond_begin','bond_end','order','stereo','direction','end_role']);for(const bond of view.bonds)rows.push([bond.begin,bond.end,bond.order,bond.stereo??'',bond.direction??'',bond.end_role??'']);rows.push([]);}for(const plot of doc.results.plots||[]){rows.push(['plot',plot.title]);for(const s of plot.series||[]){rows.push([s.name,plot.x_label,plot.y_label]);s.x.forEach((x,i)=>rows.push(['',x,s.y[i]]));}rows.push([]);}for(const table of doc.results.tables||[]){rows.push(['table',table.title]);rows.push(table.columns);rows.push(...table.rows);rows.push([]);}return rows.map(r=>r.map(v=>`"${String(v??'').replaceAll('"','""')}"`).join(',')).join('\r\n');}

async function importFile(file){
  if(!file)return;if(file.size>4_000_000){notify('Choose a JSON document smaller than 4 MB.',true);return;}
  try{const document=await api('/api/exchange/validate',JSON.parse(await file.text()));
    if(document.schema==='polymer-material/1.0'){S.material=document;clearSpecimenInputs();updateSample();updateExportState();setView('material');notify('Material record imported. Recalculate to associate results with this record.');}
    else if(document.schema==='polymer-session/1.0'){S.samples=new Map((document.samples||document.runs.map(r=>r.material)).map(m=>[m.sample_id,m]));S.material=document.material;S.runs=document.runs;S.results={};S.inputs={};S.dirty={};for(const r of S.runs.filter(r=>r.material.sample_id===S.material.sample_id)){S.results[r.application.id]=r;S.inputs[r.application.id]=clone(r.inputs);S.dirty[r.application.id]=true;}selectApp(S.active);updateSample();setView('workflow');notify('Workspace restored. Imported calculation checksums were verified. Choose View result to restore the exact sample associated with a past run.');}
    else {S.material=document.material;clearSpecimenInputs();markAllDirty();S.runs.push(document);S.runs=S.runs.slice(-50);S.results[document.application.id]=document;S.inputs[document.application.id]=clone(document.inputs);S.dirty[document.application.id]=false;updateSample();if(S.apps.some(a=>a.id===document.application.id))selectApp(document.application.id);else setView('workflow');notify('Calculation imported. Use Workflow to transfer compatible inputs to another laboratory.');}
  }catch(e){notify(`Import failed: ${e.message}`,true);}finally{$('import-file').value='';}
}

async function boot(){
  document.querySelectorAll('[data-view]').forEach(b=>b.addEventListener('click',()=>setView(b.dataset.view)));
  $('material-nav').addEventListener('click',()=>setView('material'));$('workflow-nav').addEventListener('click',()=>setView('workflow'));
  $('model-form').addEventListener('submit',e=>{e.preventDefault();run();});
  $('preset').addEventListener('change',()=>{if($('preset').value==='')return;const preset=current().presets[Number($('preset').value)];S.inputs[S.active]=Object.fromEntries(current().fields.map(f=>[f.key,clone(f.default??null)]));Object.assign(S.inputs[S.active],preset.inputs);S.dirty[S.active]=true;$('preset-description').textContent=preset.description;renderFields();updateExportState();$('run-state').textContent='Scenario loaded · review and run';});
  $('export-result').addEventListener('click',()=>{const d=S.results[S.active];exportJSON(`${slug(d.material.sample_id)}-${d.application.id}.json`,d);});
  $('export-csv').addEventListener('click',()=>{const d=S.results[S.active];download(`${slug(d.material.sample_id)}-${d.application.id}.csv`,csvFor(d),'text/csv');});
  $('print-btn').addEventListener('click',()=>window.print());
  $('session-btn').addEventListener('click',async()=>{try{const session=await api('/api/exchange/validate',{schema:'polymer-session/1.0',material:S.material,runs:S.runs,samples:Array.from(S.samples.values())});exportJSON('polymer-workspace.json',session);}catch(e){notify(`Workspace export failed: ${e.message}`,true);}});
  $('import-btn').addEventListener('click',()=>$('import-file').click());$('import-file').addEventListener('change',()=>importFile($('import-file').files[0]));
  try{const catalog=await api('/api/catalog');S.apps=catalog.apps;S.template=clone(catalog.material);S.material=catalog.material;updateSample();
    if(!S.apps.length){$('app-description').textContent='Install at least one laboratory package in this environment to begin.';notify('No laboratory packages are installed. See the README installation instructions.');return;}
    for(const app of S.apps){const b=el('button','side-link');b.dataset.app=app.id;b.append(el('span','nav-symbol',symbols[app.id]));const text=el('span',null,shortNames[app.id]);text.append(el('small',null,app.id==='chain'?'Population · conformation':app.id==='transport'?'Sorption · diffusion':app.id==='network'?'Equilibrium · release':app.id==='chemistry'?'Repeat units · atom bonds':app.id==='synthesis'?'Stoichiometry · sequence':'Correlations · conductivity'));b.append(text);b.addEventListener('click',()=>selectApp(app.id,true));$('app-nav').append(b);}
    selectApp(catalog.selected||S.apps[0].id,true);
  }catch(e){notify(e.message,true);$('app-description').textContent='The local workspace could not be initialized.';}
}
boot();
