'use strict';
// Apply an already-frozen PS/cyclohexane model to new studies. This script never fits a coefficient.
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');
const digest = data => crypto.createHash('sha256').update(data).digest('hex');
const frozenDigest = 'c7d00b3e182755bca3fa6c1c2b8032a46060d2e30a386d6bc32b577fbe5e626d';
const targetCondition = {chemistry:'PS', family:'vinylaromatic', solvent:'cyclohexane', T_C:34.5, mass_basis:'Mw'};
const localFormat = value => value.toLocaleString('tr-TR',{maximumFractionDigits:2});
const avg = vals => vals.reduce((acc,v)=>acc+v,0)/vals.length;
function checkFreeze(frozenPath, inputPath, original) {
  const raw = fs.readFileSync(frozenPath);
  if (digest(raw)!==frozenDigest) throw new Error('Dondurulmuş model sürümü veya dosya bütünlüğü değişti.');
  const archive=JSON.parse(raw);
  if(digest(fs.readFileSync(inputPath))!==archive.oldDataSHA256) throw new Error('56 eski kaydın veri sürümü değişti.');
  const reference=original.find(r=>r.id===archive.reference.id);
  if (!reference || Object.entries(archive.reference).some(([key,val])=>reference[key]!==val)) throw new Error('Dondurulmuş referans değişti.');
  return archive;
}
function freezePredict(mass, frozen, x) {
  const m=frozen.model, ref=frozen.reference;
  const logAt=m0=>m.fixed.b + 0.5*Math.log(m0/100000) + 0.5*m.fixed.slope*(x-m.fixed.center);
  const baseline=Math.exp(logAt(mass)), control=ref.eta_mL_g*Math.sqrt(mass/ref.mass_g_mol);
  const adaptive=Math.exp(logAt(mass)+m.agreement.retention*(Math.log(ref.eta_mL_g)-logAt(ref.mass_g_mol)));
  return {baseline,control,adaptive};
}
function inspect(row, source, oldRows) {
  const flags=[];
  if (row.chemistry!==targetCondition.chemistry || row.family!==targetCondition.family || row.solvent!==targetCondition.solvent) flags.push('wrong_chemistry_or_solvent');
  if (![row.T_C,row.mass_g_mol,row.eta_mL_g].every(Number.isFinite) || row.mass_g_mol<=0 || row.eta_mL_g<=0) flags.push('invalid_measurement');
  if (!source || !row.doi || !row.id || row.doi!==source.doi || row.mass_basis!==source.mass_basis || row.T_C!==source.T_C) flags.push('invalid_provenance');
  if (oldRows.some(r=>r.id===row.id || r.source_doi===row.doi)) flags.push('overlaps_old_row_or_doi');
  if (source?.knownOverlap || row.knownOverlap) flags.push('documented_sample_reuse');
  if (row.T_C!==targetCondition.T_C) flags.push('temperature');
  if(row.mass_basis!==targetCondition.mass_basis) flags.push('mass_basis');
  if(flags.length===0 && !(source?.strictEligible && source?.independence_note && source?.primary_location)) flags.push('unverified_independent_origin');
  return flags;
}
function run(test, oldRows, frozen, descriptor) {
  if (test.units?.eta_mL_g!=='mL/g' || test.units?.mass_g_mol!=='g/mol' || test.units?.T_C!=='°C') throw new Error('Dış kaynakta açıkça eşlenmiş birimler gerekli.');
  if (new Set(test.rows.map(r=>r.id)).size!==test.rows.length) throw new Error('Dış ölçüm kimliği yinelenmiş.');
  const nativePS=oldRows.filter(r=>r.chemistry==='PS');
  const train=oldRows.filter(r=>frozen.model.trainedIds.includes(r.id));
  const interval=rows=>[Math.min(...rows.map(r=>r.mass_g_mol)),Math.max(...rows.map(r=>r.mass_g_mol))];
  const modelMassRange=interval(train),sourcePSRange=interval(nativePS);
  const inRange=(mass,range)=>mass>=range[0]&&mass<=range[1];
  const out={frozen_sha256:frozenDigest,rule:'Yeni katsayı uyumu yapılmadı. Ölçümler kimya, çözücü, nominal sıcaklık, mol kütlesi tanımı ve kaynak kökeniyle tarandı. Uyuşmazlıkları olan ayrı duyarlılık karşılaştırmaları ana doğrulama paydasından çıkarılır.',
    condition:targetCondition, modelMassRange,sourcePSRange, reference:frozen.reference,studies:[],counts:{strict:0,sensitivity:0,excluded:0},rows:[]};
  for(const row of test.rows) {
    const src=test.sources[row.study],flags=inspect(row,src,oldRows);
    const soft=flags.length===1 && (flags[0]==='temperature' && Math.abs(row.T_C-targetCondition.T_C)<=0.30000001 || flags[0]==='mass_basis' && row.mass_basis==='M_supplier');
    const status=flags.length===0?'strict':soft?'sensitivity':'excluded';
    const record={...row,status,flags,expected:targetCondition,modelMassDomain:inRange(row.mass_g_mol,modelMassRange)?'interpolation':'extrapolation',databasePSMassDomain:inRange(row.mass_g_mol,sourcePSRange)?'interpolation':'extrapolation'};
    if(status!=='excluded'){
      // For sensitivity studies the numeric source mass is left intact. We evaluate the frozen
      // function at that number, without mapping supplier mass to Mw or correcting temperature.
      Object.assign(record,freezePredict(row.mass_g_mol,frozen,descriptor));
      record.errors=Object.fromEntries(['baseline','control','adaptive'].map(k=>[k,{signed_mL_g:record[k]-row.eta_mL_g,
        absolute_mL_g:Math.abs(record[k]-row.eta_mL_g),signedPct:100*(record[k]/row.eta_mL_g-1)}]));
    }
    out.rows.push(record);out.counts[status]++;
  }
  for(const [study,meta] of Object.entries(test.sources)) {
    for(const status of ['strict','sensitivity','excluded']) {
      const records=out.rows.filter(r=>r.study===study && r.status===status);
      if(!records.length)continue;
      const stats={study,status,label:meta.label,doi:meta.doi,n:records.length,reportedMassBasis:meta.mass_basis,reportedTC:meta.T_C,
        flags:[...new Set(records.flatMap(r=>r.flags))],sourceLocation:meta.table||meta.primary_location, sourceEvidence:meta.evidence||meta.independence_note,independenceNote:meta.independence_note||null,
        observedMassRange:interval(records),outsideTrainingMass:records.filter(r=>r.modelMassDomain==='extrapolation').length,outsideDatabasePSMass:records.filter(r=>r.databasePSMassDomain==='extrapolation').length};
      if(status!=='excluded')stats.scores=Object.fromEntries(['baseline','control','adaptive'].map(k=>[k,{
        mape:avg(records.map(r=>Math.abs(r.errors[k].signedPct))),bias_mL_g:avg(records.map(r=>r.errors[k].signed_mL_g)),mae_mL_g:avg(records.map(r=>r.errors[k].absolute_mL_g))}]));
      out.studies.push(stats);
    }
  }
  return out;
}
function normalize(results,sources) {
  return results.studies.filter(s=>s.status!=='excluded').map(study=>{
    const strict=study.status==='strict';
    const reason=study.flags.includes('temperature')?'Kaynak 34,8 °C, dondurulmuş model 34,5 °C. Sayısal fark sıcaklık etkisini de içerir.':
      study.flags.includes('mass_basis')?'Sıcaklık 34,5 °C. Kaynakta üretici mol kütlesinin ortalama türü belirsiz. Sayısal kütle dönüştürülmeden kullanıldı.':'Makaledeki koşullar ve mol kütlesi tanımı dondurulmuş modele uyuyor.';
    const performance=study.scores;
    return {id:'theta-external-'+study.study+'-'+study.status,label:(strict?'Yeni yayın · ':'Duyarlılık · ')+study.label,
      description:strict?'Bağımsız yayındaki ölçümlerle, dondurulmuş model ve tek standart karşılaştırması.':reason,
      baseline:'Dondurulmuş yapı modeli (Ridge)',candidate:'Dondurulmuş tek standartlı θ modeli',metric:'MAPE',unit:'mL/g',aggregation:'pooled',
      scope:study.n+' fraksiyon · '+study.label+' · PS/siklohekzan · '+localFormat(study.reportedTC)+' °C · '+study.reportedMassBasis,
      method:'Baysal PS2 referansı (90 kg/mol, 24,8 mL/g) ve eğitilmiş düzeltme katsayısı bağımsız kaynaklardan önce donduruldu. Kaynağın ham mol kütlesi ve [η] tablolardan alındı.',
      limit:(!strict?'Ana bağımsız doğrulama sayısına girmeyen duyarlılık hesabı. ':study.independenceNote+' ')+reason+' Öğrenilmiş kısmi düzeltmenin bu satırlardaki MAPE değeri %'+localFormat(performance.adaptive.mape)+'.',
      links:[{label:study.label,url:sources[study.study].url},{label:'PS2 referansı · Baysal ve Uyanık',url:'https://doi.org/10.1016/0032-3861(92)90695-S'}],
      rows:results.rows.filter(r=>r.study===study.study && r.status===study.status).map(r=>({id:r.id,group:study.label,label:r.sample,observed:r.eta_mL_g,
        baseline:r.baseline,candidate:r.control,adaptive:r.adaptive,physical:r.control,
        condition:'Kaynak '+localFormat(r.T_C)+' °C · '+r.mass_basis+' '+localFormat(r.mass_g_mol/1000)+' kg/mol · model 34,5 °C, M𝑤',
        sourceUrl:sources[study.study].url,sourceLabel:study.label,referenceId:results.reference.id,referenceUrl:'https://doi.org/'+results.reference.source_doi}))};
  });
}
if(require.main===module){
  const root=path.resolve(__dirname,'..');
  const web=process.argv[2]?path.resolve(process.argv[2]):fs.existsSync(path.join(root,'dist/family-engine.js'))?path.join(root,'dist'):root;
  const inputDir=path.join(web,'family-inputs');
  const input=process.argv[3]||path.join(inputDir,'external-theta-ps.json');
  const manifest=process.argv[4]||path.join(inputDir,'theta-ps-frozen.json');
  const destination=process.argv[5]||inputDir;
  const rawFile=path.join(web,'family-intrinsic-data.js');
  const base=require(rawFile).rows, frozen=checkFreeze(manifest,rawFile,base);
  const data=JSON.parse(fs.readFileSync(input,'utf8'));
  const results=run(data,base,frozen,require(path.join(web,'family-engine.js')).descriptor('PS').logBondsPerMass);
  results.external_data_sha256=digest(fs.readFileSync(input));
  const dump=(name,doc)=>fs.writeFileSync(path.join(destination,name),JSON.stringify(doc,null,2)+'\n');
  fs.mkdirSync(destination,{recursive:true});
  dump('external-theta-screening.json',results);
  dump('external-theta-protocols.json',normalize(results,data.sources));
  console.log(JSON.stringify({counts:results.counts,studies:results.studies.map(s=>({label:s.label,n:s.n,status:s.status,scores:s.scores}))}));
}
module.exports={inspect,run,normalize,checkFreeze};
