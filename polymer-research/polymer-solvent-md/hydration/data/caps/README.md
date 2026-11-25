# PEtOx trimerinin uç grupları

Modelde ilk amid azotuna bir metil, son omurga karbonuna bir hidroksil bağlıdır. `CN(C(=O)CC)CCN(C(=O)CC)CCN(C(=O)CC)CCO` kimliği α ucundan ω ucuna doğru yazılmış sonlu moleküldür. EtOx'un metil tosilatla başlatılması ve NaOH ile sonlandırılmasıyla α-metil/ω-hidroksil PEtOx hazırlanması deneysel olarak bildirilmiştir. Simülasyonda bu uç gruplarına sahip nötr, tam üç tekrarlı tür seçilmiştir. [Deneysel uç grup çalışması](https://doi.org/10.1039/D2PY00437B)

| Model | Formül | α ucundaki amid | ω ucu | Kuvvet alanı kütlesi (Da) | 500 su ile oligomer kütle kesri |
|---|---|---|---|---:|---:|
| H/H | C15H29N3O3 | Sekonder N–H | –CH₂CH₃ | 299.412 | 0.0321704 |
| Me/OH | C16H31N3O4 | Tersiyer N–CH₃ | –CH₂CH₂OH | 329.438 | 0.0352827 |

`source/chemical_model.json` polimerin açık kimliğini ve özgün Amber topolojisindeki atom sırasını içerir. `source/source_model.json` konformer taramasını ve seçilen başlangıç geometrisini saklar. `parameters/parameterization.json` kaydı kuvvet alanı, yazılım sürümleri, yük dengelemesi, dosya özetleri ve atom tipleme denetimlerini içerir. Koordinatlar, kısmi atom yükleri ve topoloji aynı dizindedir.

Me/OH zinciri 54 açık atom ve 53 kovalent bağ içerir. Kanonik atom sırasında α, orta ve ω tekrarlarına ait karbonil oksijenlerinin sıfır tabanlı indisleri sırasıyla 3, 20 ve 16'dır. ω-alkol oksijeni 14, hidrojeni 43'tür. Bu atomlarla N–metil karbonu, karbonil C, amid N ve son omurga karbonu `chemical_model.json` içinde ayrı ayrı belirtilmiştir. Bağ grafı, açık hidrojenler ve atom sırası Amber topolojisiyle denetlenmiştir.

Her iki koşulda üç başlangıç, aynı su sayısı, sıcaklık, basınç, kuvvet alanı ailesi, gözlem tanımları ve zaman protokolü kullanılmıştır. İki uç birlikte değiştiği için tek bir karşılaştırma uçların katkısını ayrı ayrı tayin edemez. Hidroksil ucunda hem suya verilen hem sudan alınan hidrojen bağları ölçülebilir. Toplam kütle ağırlıklı Rg ayrıca uçların kütlelerinden etkilenir. Bu nedenle ortak dokuz omurga ağır atomunun Rg'si ve α-amid N ile ω-omurga C arasındaki uzaklık da hesaplanır.

Proje ana klasöründen yeni bir vaka oluşturmak için aşağıdaki komutlar kullanılabilir. Üç tohumun her biri için ilk MD tamamlandıktan sonra son konfigürasyondan uzatma başlatılır.

```bash
python hydration/experiments/capped_petox.py --base . --out caps_case
for seed in 13031 27011 49033; do
    python hydration/experiments/capped_petox.py --base . --out caps_case --seed "$seed"
    python hydration/experiments/extend_hydration.py --base . --params caps_case/parameters --source caps_case/pilot --out caps_case/run --seed "$seed" --discard-ps 12 --production-ps 120 --stride-ps 1
done
```
