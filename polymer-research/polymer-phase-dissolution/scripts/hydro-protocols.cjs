'use strict';
const fs=require('node:fs'),path=require('node:path');
function present(data) {
  const names={pbd:'1,4-Polibütadien',pi:'1,4-Poliizopren',pib:'Poliizobütilen'};
  const lit=[{label:'Fetters ve ark. · deney tabloları 2, 4 ve 6',url:'https://doi.org/10.1063/1.555949'}];
  const stackLit=[...lit,{label:'Breiman · çapraz doğrulamayla model birleştirme',url:'https://doi.org/10.1023/A:1018046112532'}];
  function common(r,main,alternative,context){
    return {id:r.id,group:names[r.family],label:r.id,observed:r.eta_mL_g,baseline:r[main],candidate:r[alternative],
      condition:context,sourceLabel:r.source_doi,sourceUrl:'https://doi.org/'+r.source_doi};
  }
  function anchorRow(r,main,alternative) {
    const locale=num=>num.toLocaleString('tr-TR',{maximumFractionDigits:3});
    return common(r,main,alternative,'Siklohekzan · '+locale(r.M_g_mol/1000)+' kg/mol · referans '+r.anchor_id+
      ' ('+locale(r.anchor_M/1000)+' kg/mol, '+locale(r.anchor_eta)+' mL/g, DOI '+r.anchor_doi+')');
  }
  const referenceRows=data.references.groups.flatMap(g=>g.rows);
  const baseReference={metric:'MAPE',unit:'mL/g',aggregation:'macro',
    scope:'Aynı kimya ve çözücü için yayımlanmış en az bir [η] ölçümü gerekir. Hedef ölçüm modele verilmez.',
    method:'Mol kütlesi üssü yalnız diğer polimer kimyalarından öğrenilir. Aynı kimya ve ortamda farklı bir yayından kayıtlar arasında log mol kütlesi en yakın olan referans seçilir. Viskozite [η]ᵣₑ𝒻 (M/Mᵣₑ𝒻)ᵃ ile hesaplanır.',
    baseline:'Mol kütlesi modeli',candidate:'Tek referansla aktarım',links:lit};
  const main={...baseReference,id:'hydro-anchor',label:'Siklohekzanda bir yayımdan diğerine aktarım',
    description:'Bir polimerin yayımlanmış viskozitesi, başka bir çalışmadaki numuneyi daha iyi tahmin etmeyi sağlıyor mu?',
    limit:'18 hedef numune, iki kimya ve dört özgün yayın vardır. Aynı DOI kaynak ve hedefte kullanılamaz. Referans kaydı farklı hedeflerde yeniden kullanılabilir. Derleme, yakın deney sıcaklıklarını nominal 25 °C altında birleştirir. M türü poliizoprende derleme tanımıyla korunmuştur.',
    rows:referenceRows.map(r=>anchorRow(r,'mass','massCalibrated'))};
  const external={...baseReference,id:'hydro-colby',label:'Polibütadien · yeni bir yayındaki sekiz numune',description:'Sekiz polibütadien ölçümü, literatürden bir referansla aktarılan modele ayrı bir çalışmada sınama sağlıyor.',
    limit:'Bir yayın, bir yeni hedef serisi ve bilinen bir polimer. Hedef koşul siklohekzanda 25 °C’dir. Dört kayıtlı PBD referansından kütleye göre iki farklı kayıt kullanıldı. Referanslar ve hedeflerin sayısal girdi aralıkları farklıdır.',
    links:[...lit,{label:'Colby ve ark. · hedef [η], Tablo II',url:'https://doi.org/10.1021/ma00175a030'},{label:'Roovers ve Martin · referans numuneler',url:'https://doi.org/10.1002/polb.1989.090271209'}],
    rows:data.externalStudy.rows.map(r=>anchorRow(r,'mass','massCalibrated'))};
  const eb={id:'hydro-size-prior',label:'Zincir boyutunda Bayesçi düzeltme',description:'Aynı numunelerde, bir referansın hidrodinamik katsayısı fiziksel modele ne kadar katkıda bulunuyor?',metric:'MAPE',unit:'mL/g',aggregation:'macro',baseline:'Hacim modeli (Rg³/M)',candidate:'Bayesçi tek referans',
    scope:'18 numune, iki kimya, siklohekzan. Hedef ve referans için ölçülmüş Rg gereklidir.',
    method:'ln([η] M/Rg³) düzeltmesi, eğitimdeki diğer kimyalardan öğrenilen grup içi ve gruplar arası varyansla küçültülür. Yakın referansın DOI’si hedef DOI’sinden farklıdır.',
    limit:'Her dış katta varyans yalnız iki eğitim kimyasından kestirilir. Poliizopren kazancı poliizobütilende tekrarlanmaz.',
    links:lit,rows:referenceRows.map(r=>anchorRow(r,'size','sizeBayes'))};
  const negative={id:'hydro-stacking',label:'Üç yeni kimyada model birleştirme',description:'Mol kütlesi ve ölçülen zincir boyutuna dayalı iki modelin bileşimi görülmemiş bir kimyaya aktarılıyor mu?',metric:'MAPE',unit:'mL/g',aggregation:'macro',baseline:'Mol kütlesi modeli',candidate:'İç katlarda öğrenilen bileşim',
    scope:'22 ölçüm, üç ayrı kimya, siklohekzan. Her turda bir kimyanın tüm numuneleri dışarıda.',
    method:'İki modelin log tahminleri konveks ağırlıkla birleştirilir. Ağırlık dış eğitim kümesinde birer kimya ayrılarak öğrenilir. Dış hedeflerin [η] değerleri ağırlığa girmez.',
    limit:'Öğrenilen ağırlık poliizopren için 1, poliizobütilen için 0 oldu. Az kimyayla ağırlık kararsızdır. Referanssız temel model bu sınamada daha düşük ortalama hata verir.',
    links:stackLit,rows:data.folds.filter(f=>f.split==='chemistry').flatMap(f=>f.rows.map(r=>common(r,'mass','stack','Siklohekzan · nominal 25 °C · Rg '+r.Rg_nm+' nm · boyut ağırlığı '+r.sizeWeight.toFixed(3))))};
  return [main,external,eb,negative];
}
if(require.main===module){
  const input=path.resolve(process.argv[2]||'dist/hydro-stack-results.js');
  fs.writeFileSync(process.argv[3]||'protocols.json',JSON.stringify(present(require(input)),null,2)+'\n');
}
module.exports={present};
