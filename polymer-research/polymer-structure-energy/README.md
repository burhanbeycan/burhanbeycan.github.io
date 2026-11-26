# Polymer Structure & Energy

Bu bağımsız uygulama, kullanıcı tarafından tanımlanan kimyasal yapıdan sonlu bir zincir oluşturur, bağlanma düzenini ve stereokimyayı denetler, izole konformerler için klasik potansiyel enerji hesaplar. Ekrandaki sayılar kendi koşularımızdan gelir. Yapı birleştirme, taktisite denetimi, kayıt/önbellek mantığı, uygulama ve analiz akışı bu projede yazılmıştır. RDKit ETKDGv3 ve MMFF94s konformer ve enerji motorlarıdır. Atom editörü görselleştirme için OpenChemLib kullanır.

## Bilgisayarda çalıştırma

Python 3.12 ile proje klasöründe çalıştırın.

```bash
python -m venv .venv
```

macOS/Linux için `source .venv/bin/activate`, Windows PowerShell için `.venv\Scripts\Activate.ps1` ile ortamı etkinleştirin.

```bash
python -m pip install -r requirements.txt
python app.py
```

Tarayıcıda [yerel uygulamayı açın](http://127.0.0.1:8765). `python app.py --help` bağlantı noktası ve kayıt dosyası seçeneklerini listeler. Numune girilmeden herhangi bir örnek yüklenmez.

## Zinciri oluşturma

| Girdi | Oluşturulan kimyasal nesne |
| --- | --- |
| Tek, tam SMILES veya çizim | Açıkça verilmiş sonlu graf, tanımlı dallar ve stereokimya ile korunur |
| Tekrar birimi ve DP | İki tek bağla bağlı `[*:1]` ve `[*:2]` üzerinden baş–kuyruk doğrusal homooligomer. İki serbest uç H ile kapatılır |
| A/B birimleri ve açık dizi | İki farklı polimerleşmiş birimin girilen sırayla bağlı olduğu tek zincir |
| Yıldız çekirdek ve kol birimi | 3–6 numaralı uca sahip çekirdeğe eş uzunlukta kollar bağlanır. Dış uçlar H ile kapatılır |

Bileşim, zincir mimarisi ve taktisite bağımsız olarak kayda girer. Tekrar biriminin üzerindeki yan grup, tek başına makromolekülün dallı olduğunu göstermez. Birim girişleri polimerleşmiş tekrar birimleridir. Otomatik izotaktik/sindiyotaktik/ataktik seçim, bir adet uygun –CH₂–CH(R)– bölgesi olan, aynı yönde birleşen vinil tekrarlarını kapsar. Sindiyotaktik seçimde iç birimlerin yerel bağ yönüne göre *racemo* ikilileri doğrulanır. H ile kapatmanın kaldırdığı stereomerkezler sayılmaz. Ataktik seçimi belirli tohumla üretilmiş tek bir dizi gerçekleştirir. Bu dizi numunenin taktisite dağılımı değildir.

İlk sürümde doğrusal zincir veya bir kol için DP 1–32, kimlik için üst sınır 150 ağır atom, enerji hesabında üst sınır 100 ağır atomdur. Uç gruplarını, düzensiz dalları veya başka stereokimya örüntülerini tam sonlu SMILES ile tanımlayabilirsiniz. Sonsuz ağlar, numunenin molar kütle dağılımı ve sentezlenebilirlik ayrı çalışma gerektirir.

## Hesap kapsamı

Varsayılan tarama sabit tohumla 12 ETKDGv3 adayı üretir, açık H ile MMFF94s minimizasyonu yapar. Yapı, sonlu koordinat, stereokimya ve yakınsama denetimlerini geçen adaylar raporlanır. Her adayın ΔE değeri aynı yapıdaki kabul edilen örneklerin en düşük enerjisine göredir. Farklı adaylar aynı minimuma ulaşabilir. Sıralama sıcaklıkta denge popülasyonu vermez. Model enerjisinden sentez kolaylığı, izotaktik–sindiyotaktik tercih veya çözücü uyumu çıkarılmaz. İzole koordinatlardan hesaplanan kütle ağırlıklı Rg bu tek konformere aittir.

## Sonraki uygulamaya aktarım

| Kayıt | Tekrar kullanım koşulu |
| --- | --- |
| `recipe_id` | Yapım yolu, girdi birimleri, sıra, DP/kol sayısı ve varsa stereo tohumu |
| `structure_id` | Aynı kanonik izomerik sonlu graf. Başka bir yapım yoluyla aynı moleküle ulaşılabilir |
| `energy.model_id` | Aynı graf, yazılım sürümü ve hesap protokolü. Farklı tohum veya yöntem yeni kayıttır |
| `artifact_id` | Bu yapım kaydının ve hesaplanan sonucun birlikte içerik kimliği |

`POST /api/structure` kimyasal kimliği, `POST /api/analyze` yapı ve enerji sonucunu döndürür. İkisi de aynı JSON girdi biçimini kullanır. Sonuç ekranından JSON indirilebilir. `GET /api/records/<artifact_id>` kalıcı kayıt kopyasını getirir. Geçici `cache_hit` ve `cache_source` alanları o talebe aittir. Varsayılan kayıt yeri `~/.polymer-research-studio/structure-energy.sqlite3` dosyasıdır.

Örneğin aşağıdaki girdiyi `zincir.json` olarak kaydederek aynı işlemi komut satırında çalıştırabilirsiniz.

```json
{"mode":"repeat","smiles":"[*:1]CCN(C(=O)CC)[*:2]","dp":3,"tacticity":"as_drawn"}
```

```bash
python app.py --spec zincir.json
python -m unittest discover -s . -p "test_*.py"
```

## Çözücü bileşenleri

Su, etanol, asetonitril, aseton, THF, DMF ve DMSO kendi ETKDGv3/MMFF94s akışımızla taranmıştır. `GET /api/solvents` tüm kayıtları, `GET /api/solvents/water` tek kaydı verir. `python solvent_registry.py --write solvent_components.json` sayıları yeniden üretir. Aynı kimyasal graf ve protokol istendiğinde dağıtılan sonuç doğrudan kullanılır. Yazılım sürümü veya ayarlar farklıysa yeni hesap açılır. Kaydın 3B koordinatları başlangıç geometrisi sunar. MD topolojisi, sıvı kutusu, su modeli ve yük modeli kullanılan kuvvet alanına uygun ayrıca tanımlanır.

## Saf sıvı kaydı

Çözücü sekmesindeki saf su pilotu, tekil H₂O konformerinden bağımsız bir OpenMM hesabıdır. Klasik TIP3P modeliyle periyodik 216 moleküllük kutu, PME, 298,15 K ve 1 bar kullanılır. Üç farklı rastgele tohumla başlatılan her koşuda 12 ps dengelemeden sonra 32 ps üretim yapılır. Gösterilen yoğunluk, üretim bölümünde kutu kütlesi ve anlık hacimden hesaplanır. Koşu ortalamalarının saçılımı ve blok ortalamaları, kısa örneklemenin değerlendirilmesi için saklanır. Üretim çerçevelerinden periyodik kutunun potansiyel enerjisi su sayısına bölünerek mol su başına ortalama potansiyel enerji (kJ/mol) de verilir. Enerji bu modelin referansına ve protokolüne özgüdür. Bu pilot saf sıvıya aittir.

`GET /api/liquid-reference` sonuçları ve koşu kimliğini sunar. `pure_liquid/` dizinindeki betik, girdi ayarları, zaman serisi ve son periyodik durumlar bağımsız olarak kullanılabilir. Su modeli veya termodinamik koşullar değişince yeni bir protokol gerekir. Çözünen eklendiğinde karışımın ayrıca dengelenmesi gerekir. Bu uygulamanın tüm yöntemleri klasik fizik ve kimyasal graf hesaplarıdır.

## Önerilen bağımlılık sırası

1. Bu uygulama gerçek kimyasal grafı, model zincirini ve başlangıç konformerini sağlar. Saf çözücü kütüphanesi bu işten bağımsız hazırlanabilir.
2. Bir sonraki polimer–çözücü uygulaması iki bileşenin uygun klasik parametrelerini seçer. Nihai sistemi seçilen bileşim, sıcaklık, basınç ve fazda kurar ve dengeler. Tekrarlar, sonlu-boyut etkisi ve örnekleme kontrol edilir. Aynı koşu yörüngesi uyumlu Rg dağılımı, radyal dağılım ve hareket analizlerine girdi olabilir.
3. Serbest enerji gerekiyorsa hedef süreç ve referans hâli seçilir, ona özgü örnekleme yapılır. Kuşkulu, eksik veya uyumsuz parametre durumunda sonuç alanı gerekçeyle boş kalır.
4. AI için uygun ilk görev, doğrulanmış kendi koşularımızın etiketlediği veriyle benzer adayları sıralamak ve yeni MD koşusu seçmektir. Bağımsız polimer kimyalarında hata ve uygulama alanı belirlendikten sonra, ML kestirimleri açıkça ayrı bir yöntem olarak raporlanabilir.

### Yöntem kaynakları

[RDKit konformer üretimi](https://www.rdkit.org/docs/RDKit_Book.html) · [RDKit MMFF kuvvet alanı](https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html) · [IUPAC sindiyotaktik makromolekül](https://goldbook.iupac.org/terms/view/S06224) · [GROMACS solvasyon serbest enerjisi](https://tutorials.gromacs.org/docs/awh-free-energy-of-solvation.html) · [OpenMM simülasyon yöntemleri](https://docs.openmm.org/latest/userguide/application/02_running_sims.html)

Paketlenmiş OpenChemLib JS 9.26.0 için BSD-3-Clause lisansı `assets/ocl-vendor.js` içinde yer alır.
