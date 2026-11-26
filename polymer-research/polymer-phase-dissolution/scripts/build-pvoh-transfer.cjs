// From the repository root: node scripts/build-pvoh-transfer.cjs [data-directory] [output-directory]
const fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'..');
const dir=[process.argv[2],path.join(root,'dist'),root].find(d=>d&&fs.existsSync(path.resolve(d,'pvoh-transfer-engine.js')));
if(!dir)throw Error('Cannot find the PVOH data and engine.');
const base=path.resolve(dir),out=path.resolve(process.argv[3]||path.join(base,'family-inputs'));
fs.mkdirSync(out,{recursive:true});
const E=require(path.join(base,'pvoh-transfer-engine.js')),D=require(path.join(base,'pvoh-transfer-data.js')),R=require(path.join(base,'family-pvoh.js'));
const result=E.evaluate(R.rows,D);
const report=[{
  id:'pvoh-solvent-transfer',label:'PVOH · saklı numune ve çözücüye özgü eğimler',
  description:'Her numunenin su ve DMSO ölçümleri birlikte saklanır. İki eğimin çözücüye göre değişmesine izin veren modelin hatasını eski ortak hidroliz modeli ile karşılaştır.',
  baseline:'Ortak hidroliz etkisi',candidate:'Çözücüye özgü iki eğim · kısmi havuzlama',metric:'RMSE',unit:'mL/g',aggregation:'pooled',
  scope:'25 °C · Parisi ve ark. · 5 bağımsız numune · 10 ölçüm',
  method:'Mark–Houwink kütle üssü ve hidroliz katsayısı her çözücüde uyarlanır, ardından her bir eğim ortak katsayıya doğru yarı yarıya çekilir. Değerlendirme aynı numuneye ait iki ölçümü birlikte dışarıda tutar.',
  limit:'Bu ikinci aday birincil DH denemesinden sonra incelendi. Birleşik RMSE düşse de suya ait eski model daha iyi. DMSO için eski en iyi kütle modeli 13,41 mL/g verir. Bağımsız sulu seriler ayrı kapsamda gösterilir.',
  links:[{label:'Çözücüler ve bağımsız deneyler',url:'pvoh-transfer.html'},{label:'Parisi · birincil deney',url:D.urls.parisi}],
  rows:result.folds.map(r=>({id:r.grade+'/'+r.solvent,group:r.solvent==='water'?'Su':'DMSO',label:r.grade,
    observed:r.measured,baseline:r.sharedDH,candidate:r.partialBoth,
    condition:'M𝑤 = '+r.Mw_g_mol+' g/mol · DH %'+r.DH_pct+' · 25 °C',sourceUrl:D.urls.parisi,sourceLabel:'Parisi, Tablo III'}))
}];
const serialize=(name,x)=>fs.writeFileSync(path.join(out,name),JSON.stringify(x,null,2)+'\n');
serialize('pvoh-transfer-results.json',result);
serialize('pvoh-transfer-protocols.json',report);
