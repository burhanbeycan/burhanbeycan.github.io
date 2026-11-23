'use strict';

// Loaded on demand so that past material/result files retain their original flow.
async function renderPilots(){
  const root=$('view-pilots');root.replaceChildren();
  const head=el('div','card');head.append(el('h2',null,'LDI diurethane and hydrolysed PEtOx'),
    el('p',null,'Local dry-cluster dynamics at 298.15 K. One mixed-arm LDI diurethane was paired with a neutral EtOx dimer or its singly hydrolysed analogue. Both have a methyl α end and a hydroxyl ω end.'));
  const img=el('img','pilot-structure');img.src='/static/ldi-structure.svg';img.alt='Mixed-arm lysine diurethane with ether and ester groups. Two urethane NH donors are highlighted.';head.append(img);root.append(head);
  try{
    const data=await api('/static/ldi-blend-pilot.json');
    const controls=el('form','material-form');controls.id='pilot-results';
    const pose=recordField(controls,'pose','Starting pose',{'undirected':'Undirected approach','carbonyl contact':'Seeded N–H···O=C contact'});
    const partners=Object.fromEntries(Object.keys(data.partner_smiles).map(x=>[x,x]));
    const partner=recordField(controls,'partner','Partner sequence',partners);
    const groups=recordField(controls,'site','Acceptor site',{});
    const box=el('section','card');box.append(el('h3',null,'Contacts across independent initializations'),controls);
    const body=el('div');box.append(body);root.append(box);
    function refreshSites(){
      const run=data.runs.find(r=>r.partner===partner.value);
      const previous=groups.value;
      groups.replaceChildren();for(const c of run.contacts){const opt=el('option',null,c.group);opt.value=c.group;groups.append(opt);}
      groups.value=run.contacts.some(c=>c.group===previous)?previous:'EtOx amide carbonyl';
    }
    function draw(){
      const selected=data.runs.filter(r=>r.initialization===pose.value && r.partner===partner.value);
      const rows=selected.map(r=>{const c=r.contacts.find(x=>x.group===groups.value);return [r.seed,c.contact_count,c.frames,c.contacts_per_donor,c.occupancy_per_eligible_pair];});
      body.replaceChildren(tableCard({title:groups.value,columns:['Packing seed','Contacts summed over frames','Saved frames','Mean per host NH','Mean per eligible pair'],rows}));
      body.append(el('p',null,'Each packing is one replicate. A recorded contact satisfies N···A ≤ 0.35 nm and N–H···A ≥ 150°. Frames within a replicate are correlated.'));
    }
    pose.onchange=draw;partner.onchange=()=>{refreshSites();draw();};groups.onchange=draw;
    refreshSites();draw();
    const details=el('details','card');details.append(el('summary',null,'Method, files and interpretation'),el('p',null,data.method),
      el('p',null,'Per packing: 0.5 ps initial relaxation and 3 ps production. The two fixed-charge neutral dimers and truncated PEG/PCL arms are small-molecule controls. Equilibrium affinity and film properties require other evidence.'));
    const raw=el('a',null,'Download run-by-run results (JSON)');raw.href='/static/ldi-blend-pilot.json';raw.download='ldi-blend-pilot.json';details.append(raw);root.append(details);
  }catch(e){root.append(el('p','notice error',e.message));}
}
