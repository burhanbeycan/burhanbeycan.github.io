# Dondurulmuş PS modeli için yeni yayınların taranması

Dondurulmuş koşul lineer PS, saf siklohekzan, 34,5 °C ve ağırlık ortalamalı mol kütlesi M𝑤'dir. [η] mL/g olarak karşılaştırılır. Kaynaklardan alınan ölçümler modelin eğitiminde veya seçiminde kullanılmamıştır. Kaynak kimliği, numune kimliği, yöntem ve tekrar yayımlanan ölçümler denetlenir. Dondurulan model, katsayılar, eğitim kayıtları ve tek PS standardının kimliği `theta-ps-frozen.json` içindedir. SHA-256 değeri c7d00b3e182755bca3fa6c1c2b8032a46060d2e30a386d6bc32b577fbe5e626d'dir.

Miyaki, Einaga, Fujita ve Fukuda'nın ayrı makalesinde aynı polistiren standardından yeniden ayrılmış iki dar fraksiyon dondurulmuş koşulların tamamını sağlar. M𝑤 6.480 ve 6.550 kg/mol, deneysel [η] sırasıyla 228 ve 230 mL/g'dir. Saklı ölçümlerde MAPE, önceki model için %8,86, Baysal standardıyla fiziksel düzeltmede %7,86 ve öğrenilmiş kısmi düzeltmede %7,74'tür. İki hedef de orijinal modelin aromatik dışı eğitim kütle aralığının üzerindedir. Aynı stoktan üretilen bu iki fraksiyon tek bir deney serisidir. Üç yazar önceki Einaga çalışmasıyla ortaktır. Referans standardı başka bir gruptan gelir. Bu sınama ikinci bir laboratuvar tekrarı oluşturmaz. [Miyaki ve ark., özgün çalışma](https://doi.org/10.1021/ma60075a021) · [Fetters ve ark., Tablo 7](https://srd.nist.gov/jpcrdreprint/1.555949.pdf)

Abdel Azim ve Huglin'in beş fraksiyonunda sıcaklık 34,5 °C'dir. Kaynak kütleyi üretici değeri olarak raporlar ve ortalama türünü bildirmez. Sayısal kütle olduğu gibi kullanılarak yapılan duyarlılık karşılaştırmasında temel model MAPE'si %13,62, bir standartlı fizik kuralında %3,83 ve kısmi düzeltmede %3,71'dir. Üretici değerinden M𝑤'ye bir dönüşüm yapılmadı. 1983 ve 1985 raporları aynı seri olduğundan bir kez sayıldı. [Abdel Azim ve Huglin, birincil makale ve tablo](https://asat.journals.ekb.eg/article_26559_db5988156172c60aeea13b121939bc25.pdf)

Berry'nin 16 fraksiyonunda M𝑤 raporlanır. Tabloda verilen 34,8 °C, dondurulmuş modelin 34,5 °C koşulundan 0,3 °C daha yüksektir. 34,5 °C modelini sayısal M𝑤'de değerlendirip kaynak ölçümüyle karşılaştırınca MAPE eski yapı modeli için %15,01, bir standartlı fizik için %2,75 ve kısmi düzeltme için %2,64 olur. Sıcaklık düzeltmesi yapılmadı. Derlemelerdeki tekrar basımlar aynı gözlemlerdir. [Berry, özgün makale](https://doi.org/10.1063/1.1840854)

Bu iki duyarlılık karşılaştırması ayrı tutulur. Aralarından en küçük hatayı seçmek, doğrulamada geçerli bir dış kaynak oluşturmaz. İki seri de donmuş PS koşulunun bütün şartlarını karşılamadığından bu iki seriden kesin dış doğrulama paydasına satır eklenmez. Yalnızca koşulları tam eşleşen Miyaki serisinin iki fraksiyonu bu paydayı oluşturur. Tarama çıktısı `external-theta-screening.json` her kaynağın orijinal koşullarını, kapsam kararını, eğitim ve önceki PS veri aralıklarının dışındaki gözlemleri, numune bazında işaretli ve mutlak hataları içerir.

Yeni bir kaynağı katarken raporlanmış sayısal ölçüm, bire bir koşul, aynı M𝑤 tanımı ve bağımsız fraksiyon kökenini doğrulamak gerekir. Kütle dayanağı eksik veya farklı sıcaklıktaki duyarlılık kayıtları harici fizik tabanlı sayısal karşılaştırmada görünür. Canlı uygulama bu koşullara genişletilmez.

```sh
node scripts/build-external-theta.cjs dist
node tests/test-external-theta.cjs dist
```
