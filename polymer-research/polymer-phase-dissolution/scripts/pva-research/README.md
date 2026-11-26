# Polymer Research Studio · Ters reoloji ve zincir kimyası

Çevrimdışı arayüz için `app/pva-inverse.html` dosyasını açın. Arayüzü kullanmak için kurulum gerekmez. Polimer tasarımı, tanımlı bir ticari numunenin akışı, yaşlanma deneyi ve polyBERT ile aile-dışı sınama ayrı seçeneklerdir.

## Veriler ve deney protokolleri

`data/pva-inverse-data.csv` PVA için birincil kaynaklara bağlı birleşik çizelgedir. Her satırın kaynak çalışması, numune veya ürün kimliği, şekil/tablo konumu, özgün derişim temeli, ölçüm yöntemi ve hedefi vardır. Boş hücre bilinmeyen değerdir. Kinematik viskozite, intrinsik viskozite ve MFR bu cP hedeflerinin yerine kullanılmaz. `data/manifest.json` çalışmalara göre sayıları verir. Tek bir reometre eğrisinin onlarca noktası tek bir numuneye ait olabilir.

- `parisi/` 298,15 K’de iki çözücü için sayısallaştırılmış ölçüm işaretçilerini ve piksel doğrulamasını içerir. Suda şekil lejantı ve numune tablosu uyuşmayan 79 nokta için mol kütlesi alanı boştur. Ayrı raporlanan iki piksellik duyarlılık deneysel ölçüm belirsizliğini temsil etmez. [Parisi vd.](https://doi.org/10.1122/8.0000435)
- `model/` tek bir ticari PVA ürününün dört sıcaklık ve dört derişimde 244 kayma ölçümüne dayanır. Bu ürünün Mₙ değeri bilinmiyor. Kayma verisinden 30 s⁻¹’deki viskozite ve iki yerel akış üssü tersine kestirilir, ardından bunların derişim ve sıcaklığa bağlılığı öğrenilir. Her akış eğrisi eşit ağırlık alır. Sonuçlar numune ve laboratuvar bazında da sınanmıştır. [Carrasco-Venegas vd.](https://doi.org/10.3390/polym15071743) · [Pérez-Puyana vd.](https://doi.org/10.1088/2053-1591/aab164)
- `literature/` bağımsız sayısal ölçümler, cihaz ve hazırlama koşulları ile özgün birimlere dair kararları içerir. Dönme hızını kalibre edilmiş geometri olmadan kayma hızına veya g/20 mL değerini yoğunluk bilinmeden kütle yüzdesine çevirmeyin.
- `kinetics/` depolama yaşına bağlı viskozite ile milisaniye ölçeğindeki terminal gevşeme ölçümlerini ayrı tutar. Modelde 303,15 K ve başlangıç ölçümü bilinir. 298,00 K için zaman düzeltmesini kullanmadan önce aynı numuneye ait izotermal zaman serisi ve farklı sıcaklıklarda hız katsayısı gerekir. [Yaşlanma](https://doi.org/10.1295/koron1944.15.348) · [Terminal gevşeme](https://doi.org/10.1021/acs.macromol.5c00703)
- `descriptors/` kuru tekrar birimi, OH/ester yoğunluğu, mol kütlesi ortalamaları, koşula bağlı boyutsuz derişim ve önceden hesaplanan 600 boyutlu polyBERT vektörleri için araçları içerir. DH’den türeyen birden fazla tanımlayıcının aynı küçük veri setinde bağımsız değişken gibi davranması beklenmez. [polyBERT](https://doi.org/10.1038/s41467-023-39868-6)
- `grades/` Parisi’nin numune kimliği doğrulanmış sulu seriyle dört numune-dışı sınama yapar. Mₙ regresyonu bir karşılaştırmadır. Arayüzdeki sanal Mₙ eğrisi ise aynı DH’deki iki komşu referanstan logaritmik enterpolasyon ile oluşturulur. Bu ikincisinin bağımsız ara Mₙ doğrulaması henüz yoktur. 298,00 K düzeltmesinde başka ürünün yerel sıcaklık oranı kullanıldığından sonuç araştırma tahminidir. Numunelerin Đ değerleri de değişir.
- `transfer/` dört polimer malzemesi için kimya ve laboratuvar dışı aktarım sınamasını içerir. Sadece kararlı kayma koşullarında elde edilmiş 286 gözlem kullanılır. PVA satırları birleşik PVA çizelgesinin alt kümesidir. Derişim temelleri korunur ve aljinatın bilinmeyen kopolimer bileşimi için iki uç duyarlılık hesabı vardır. Her sıcaklık, derişim ve tuz kombinasyonundaki bir ölçüme dayalı kalibrasyon ayrı puanlanır.

## Hesabı yeniden çalıştırmak

Ana kıyaslama dosyaları Python 3, NumPy, SciPy ve pandas kullanır. Yapı denetimi için RDKit, isteğe bağlı olarak hesaplanan yeni temsiller için ilgili encoder betiğinin bağımlılıkları gerekir. Hazır vektörlerle kıyaslamada büyük model ağırlıklarını indirmeye gerek yoktur. Ham ölçüm ve sayısallaştırma kaynakları alt klasörlerde belirtilmiştir.

```sh
python data/assemble.py
python model/pva_model.py build
python model/external_product_test.py
python grades/grade_transfer.py
python kinetics/aging.py
python transfer/benchmark.py
python prepare_studio.py
node model/test-parity.cjs
```

Kayma hızına, T’ye, DH’ye veya Mₙ’ye ilişkin bilgi eksikse uygun çalışma alanını ve boş değerleri inceleyin. Yüksek doğruluk yeni bir kimya için ancak dağılım dışı veri ve uygun protokolde bağımsız reometriyle gösterilebilir.


## Sıcaklık ve kaynak aktarımı sonuçları

Aynı ticari PVA'da, her seferinde tam bir sıcaklığı dışarıda bırakan testte ters modelin eğri ağırlıklı MAPE değeri %6,04 oldu. Tam bir derişim dışarıda bırakıldığında %18,13, bağımsız Pérez-Puyana numunesi ve laboratuvarındaki desteklenen üç noktada %48,58 elde edildi. Her yeni polimer için Mₙ'ye dayalı evrensel bir doğruluk değeri mevcut değildir.

Dört PVA numunesi sırayla dışarıda bırakılan sulu veri testinde Mₙ terimli regresyonun numune ağırlıklı MAPE değeri %40,11 oldu. Hidroliz terimini eklemek bu küçük örneklemde %70,35'e çıkardı. Bu sonuç hidrolizin kimyasal olarak önemsiz olduğunu göstermez; sıcaklık, mol kütlesi, dağılım ve hidrolizin etkilerini ayırmak için dengeli bir grade dizisi gerekir. Teorik Mₙ ekranının varsayılan 60 s⁻¹ değeri, düşük ve yüksek mol kütleli %80 hidrolizli numunelerin ortak pencere bölgesindedir. Aynı ticari ürün için akış eğrisi ekranının varsayılanı 10 s⁻¹'dir.

Aile dışı aktarımda 286 ölçümden ikisi fizikokimyasal tanımlayıcılar + polyBERT modelinde ölçülen değerin ±%10'u içinde kaldı. Bu dört malzeme ve üç dış grup, laboratuvar ve protokol etkileriyle kısmen örtüşür. Yeni lineer PSMILES'i araştırma modunda kodlamak ve bağlam bilgileriyle yeni bir tahmin üretmek için `transfer/README.md` ve `descriptors/README.md` adımlarını kullanın.
