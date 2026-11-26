'use strict';
const fs = require('node:fs');
const path = require('node:path');
const root = path.resolve(__dirname, '..');
const web = process.argv[2] ? path.resolve(process.argv[2]) : fs.existsSync(path.join(root,'dist/family-engine.js')) ? path.join(root,'dist') : root;
const destination = process.argv[3] ? path.resolve(process.argv[3]) : path.join(web,'family-inputs');
const moduleDir = fs.existsSync(path.join(web,'theta-transfer-ai.js')) ? web : __dirname;
const {benchmark} = require(path.join(moduleDir,'theta-transfer-ai.js'));
const records = require(path.join(web,'family-intrinsic-data.js')).rows;
const result = benchmark(records);
fs.mkdirSync(destination,{recursive:true});
const save=(file, value)=>fs.writeFileSync(path.join(destination,file),JSON.stringify(value,null,2)+'\n');
save('theta-learning-results.json',result);
const families={polydiene:'Polidienler',polyolefin:'Poliizobütilen',vinylaromatic:'Vinil aromatikler','fetters1994-ref-46':'Einaga vd. (1979)',baysal1992:'Baysal ve Uyanık (1992)'};
const links=[{label:'Fetters vd. · veri derlemesi',url:'https://doi.org/10.1063/1.555949'},{label:'Baysal ve Uyanık · θ ölçümleri',url:'https://doi.org/10.1016/0032-3861(92)90695-S'},{label:'Einaga vd. · θ ölçümleri',url:'https://doi.org/10.1002/pol.1979.180171205'}];
const rawById=new Map(records.map(r=>[r.id,r]));
const sname={'1,4-dioxane':'1,4-dioksan','isoamyl isovalerate':'izoamil izovalerat',cyclohexane:'siklohekzan'};
const turk=n=>Number(n).toLocaleString('tr-TR',{maximumFractionDigits:2});
const common={metric:'MAPE',unit:'mL/g',aggregation:'macro',links};
function pack(base, source, comparison='candidate') {
  return {...common,...base, rows:source.points.map(p=>{
    const row=rawById.get(p.id);
    const ref=p.referenceId ? rawById.get(p.referenceId) : null;
    const condition=[sname[p.solvent],(row.temperature_basis==='approximate_compilation'?'yaklaşık ':'')+turk(p.T_C)+' °C',
      (row.mass_basis==='Mw'?'M𝑤':'Mabs')+' '+turk(p.mass/1000)+' kg/mol'];
    if(ref)condition.push('Referans '+ref.sample+' · '+turk(ref.mass_g_mol/1000)+' kg/mol · '+turk(ref.eta_mL_g)+' mL/g');
    return {id:p.id,group:families[p[source.averagingUnit]],label:p.chemistry+' · '+row.sample,observed:p.observed,
      baseline:p.baseline,candidate:p[comparison],condition:condition.join(' · '),sourceUrl:'https://doi.org/'+p.source_doi,sourceLabel:row.source, 
      ...(ref?{referenceId:ref.id,referenceUrl:'https://doi.org/'+ref.source_doi,adaptive:p.candidate,physical:p.control}:{})};
  })};
}
const protocols=[
 pack({id:'theta-one-reference',label:'θ · kayıtlı bir referans',description:'Aynı polimerin kayıtlı bir fraksiyonu, diğer mol kütlelerindeki [η] tahminini ne kadar değiştiriyor?',
   baseline:'Eski yapı modeli (Ridge)',candidate:'Tek referanslı θ ölçeklemesi',scope:'Beş hidrokarbon kimyası, üç aile, beş referans hariç 51 hedef fraksiyon.',
   method:'Her kimya için 100 kg/mol’e log ölçekte en yakın aynı koşullu fraksiyon önceden seçildi. [η] = [η]ref √(M/Mref). Hedef fraksiyonların viskoziteleri bu seçimde kullanılmadı. Ailelerin MAPE değerleri eşit ağırlık alır.',
   limit:'Standartların 51 hedefin çoğuyla aynı yayından gelmesi benzer ölçüm hatalarının aktarılmasına yol açabilir. Eğitim kimyalarından öğrenilen kısmi düzeltmede MAPE %3,10, tam referans düzeltmesinde %2,78 oldu.'},result.oneReference,'control'),
 pack({id:'theta-zero-reference',label:'θ · referanssız aile aktarımı',description:'Hedef aileden hiçbir viskozite ölçümü yokken iç doğrulama ile Ridge düzenlileştirmesini seçmek hatayı azaltıyor mu?',
   baseline:'Sabit Ridge',candidate:'İç sınamayla seçilen Ridge',scope:'56 ölçüm, üç omurga ailesi. Her dış turda bütün aile saklandı.',
   method:'Dış eğitimdeki kimyalar sırayla saklanarak beş eğim zayıflatma katsayısı karşılaştırıldı. Ağırlıklandırılmış log-hata en düşük katsayı seçildi. Aile başına ortalama mutlak yüzde hatalar eşit ağırlık aldı.',
   limit:'Tümü hidrokarbon olan beş kimya tanımlayıcı çeşitliliğini sınırlar. İç seçim eski sabit katsayının hata düzeyini iyileştirmedi.'},result.zeroShot),
 pack({id:'theta-ps-between-studies',label:'θ · iki PS yayını arasında aktarım',description:'Tek polistiren standardı başka bir çalışmadan alındığında tahmin nasıl değişiyor?',
   baseline:'Eski yapı modeli (Ridge)',candidate:'Diğer yayındaki tek standartla θ',scope:'Siklohekzan, 34,5 °C. İki bağımsız yayında toplam 12 PS fraksiyonu.',
   method:'Her hedef çalışmanın tüm PS ölçümleri referans havuzundan çıkarıldı. Diğer yayının 100 kg/mol’e log ölçekte en yakın fraksiyonu seçildi. Çalışma MAPE değerleri eşit ağırlık aldı.',
   limit:'İki kaynakta mol kütleleri farklı aralıkları örnekler. Eğitim kimyalarından öğrenilen kısmi düzeltme iki makalede ortalama %5,58, yalın fiziksel aktarım %5,62 verdi.'},result.externalPS,'control')
];
save('protocols.json',protocols);
for(const [k,v] of Object.entries(result))console.log(k,v.n,v.scores);
