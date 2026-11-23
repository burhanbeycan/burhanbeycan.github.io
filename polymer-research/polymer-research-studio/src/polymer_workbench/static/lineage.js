'use strict';

const specimenForms = {
  batch:'Synthesis batch', cast_film:'Cast film', electrospun_mat:'Electrospun mat', solution:'Solution', other:'Other'
};
const evidenceOrigins = {
  own_measurement:'Own measurement', primary_article:'Measured in the reporting article',
  secondary_citation:'Value cited from an earlier source', model_assumption:'Model assumption',
  illustrative:'Illustrative input'
};
function clearSpecimenInputs(){
  S.results={}; S.inputs={}; markAllDirty();
  if(S.active)selectApp(S.active);
}
function selectSpecimen(record){
  S.material=clone(record);
  clearSpecimenInputs();
  updateSample();
  setView('material');
}
function recordField(form,key,label,options){
  const row=el('div','field'), caption=el('label',null,label);
  const input=el(options?'select':'input');
  input.name=key;
  input.id=form.id+'-'+key;
  input.required=true;
  if(options)for(const [value,text] of Object.entries(options)){
    const item=el('option',null,text); item.value=value; input.append(item);
  }
  caption.htmlFor=input.id;
  row.append(caption,input); form.append(row); return input;
}
function specimenPanel(){
  const card=el('section','card');
  card.append(el('h3',null,'Physical specimens'));
  const table={title:'Specimen register',columns:['ID','Form','Derived from','Process'],rows:[]};
  for(const specimen of S.samples.values()){
    const lineage=specimen.lineage;
    table.rows.push([specimen.sample_id,specimenForms[lineage?.form]||'Unspecified',
      lineage?.ancestors.at(-1)||'—',lineage?.process||'—']);
  }
  card.append(tableCard(table));
  const switcher=el('form'); switcher.id='switch-specimen';
  const choice=recordField(switcher,'target','Active specimen',Object.fromEntries(
    Array.from(S.samples.values(),s=>[s.sample_id,s.sample_id+' · '+s.name])));
  choice.value=S.material.sample_id; choice.onchange=()=>selectSpecimen(S.samples.get(choice.value));
  card.append(switcher);
  const form=el('form','material-form'); form.id='new-specimen';
  const route=recordField(form,'route','Relationship',{derived:'Derived from active specimen',root:'Independent batch'});
  recordField(form,'sample_id','New specimen ID'); recordField(form,'name','Name');
  const formType=recordField(form,'form','Specimen form',specimenForms); formType.value='cast_film';
  route.onchange=()=>{formType.value=route.value==='root'?'batch':'cast_film'};
  recordField(form,'process','Preparation or processing step');
  const save=el('button','primary small','Register specimen'); form.append(save);
  form.onsubmit=async event=>{
    event.preventDefault();
    try{
      const values=Object.fromEntries(new FormData(form));
      for(const k of Object.keys(values))values[k]=values[k].trim();
      if(S.samples.has(values.sample_id))throw new Error('This specimen ID already exists.');
      if(S.samples.size>=100)throw new Error('The workspace already contains 100 specimens.');
      const {route,...node}=values;
      let child;
      if(route==='derived'){
        child=await api('/api/material/derive',{parent:S.material,...node});
      }else{
        child=clone(S.template);child.sample_id=node.sample_id;child.name=node.name;
        child.conditions={temperature_K:null,state:'unspecified',solvent:'unspecified',notes:''};
        child.chemistry.notes='';child.lineage={form:node.form,process:node.process,ancestors:[]};
        child=await api('/api/exchange/validate',child);
      }
      selectSpecimen(child);notify('Specimen registered. Enter conditions and evidence for this specimen.');
    }catch(e){notify(e.message,true);}
  };
  card.append(form);return card;
}
function evidencePanel(){
  const card=el('section','card');card.append(el('h3',null,'Evidence for '+S.material.sample_id));
  const entries=S.material.evidence||[];
  if(entries.length)card.append(tableCard({title:'Recorded evidence',columns:['Quantity','Value','Unit','Conditions / method','Origin','Reported in','Attributed to','Locator'],rows:entries.map(e=>
    [e.quantity,e.value,e.unit,e.conditions,evidenceOrigins[e.origin],e.reported_in||'—',e.attributed_to||'—',e.locator])}));
  const form=el('form','material-form');form.id='new-evidence';
  const origin=recordField(form,'origin','Evidence origin',evidenceOrigins);
  for(const [key,label] of Object.entries({quantity:'Quantity',value:'Value',unit:'Unit (1 if dimensionless)',conditions:'Method and conditions',reported_in:'Reporting article (DOI or citation)',attributed_to:'Original source cited (or unresolved)',locator:'Table, figure or laboratory record'}))recordField(form,key,label);
  function update(){
    for(const [key,show] of [['reported_in',['primary_article','secondary_citation'].includes(origin.value)],['attributed_to',origin.value==='secondary_citation']]){
      const control=form.elements.namedItem(key);control.required=show;control.disabled=!show;control.parentElement.hidden=!show;
    }
  }
  origin.onchange=update;update();
  form.append(el('button','primary small','Attach evidence'));
  form.onsubmit=async event=>{
    event.preventDefault();
    try{
      const material=clone(S.material);
      material.evidence=[...(material.evidence||[]),Object.fromEntries(new FormData(form))];
      S.material=await api('/api/exchange/validate',material);markAllDirty();updateSample();updateExportState();renderMaterial();
      notify('Evidence attached to '+S.material.sample_id+'.');
    }catch(e){notify(e.message,true);}
  };
  card.append(form);return card;
}
