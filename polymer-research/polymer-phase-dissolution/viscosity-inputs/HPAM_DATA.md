# Sulu HPAM–NaCl görünür viskozite verisi

Saeed Akbari, Syed M. Mahmood, Isa M. Tan ve Hamed Hematpour, *Journal of Petroleum Exploration and Production Technology* **8**, 887–900 (çevrimiçi 2017, cilt 2018). Dönüştürdüğümüz sayısal kayıtlar makalenin CC BY 4.0 lisansı kapsamındadır. [Yayın ve lisans beyanı](https://doi.org/10.1007/s13202-017-0375-6) · [Birincil tabloların PDF'i](https://link.springer.com/content/pdf/10.1007/s13202-017-0375-6.pdf) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Numune ve ölçüm

Ticari kısmen hidrolize poliakrilamidin (SNF Floerger) üretici tarafından bildirilen hidroliz derecesi %5 mol, mol kütlesi yaklaşık 8–9 milyon g/mol'dür. Yayında bu aralık için ortalama türü ve dispersite verilmemiştir. Toz polimer, demineralize su ve NaCl ile hazırlanmış tuzlu suda en az bir gün hidratılmıştır. Malvern Bohlin Gemini II döner reometreyle 24 °C'de artan ve azalan kayma taramalarında görünür viskozite ölçülmüştür. Çalışmada kayma histerezisi görülmediği bildirilmiştir.

## Kayıt ve dönüşüm

1. Sayfa 890'daki Tablo 1'den koşulları özgün birimlerinde aldık. Polimer derişimi ppm, NaCl derişimi g/L, kayma hızı s⁻¹ olarak saklanır. Makalenin ppm birimi kütle/hacim birimine çevrilmemiştir.
2. Sayfa 895'teki Tablo 5'in deneysel viskozite sütununu deney numarası üzerinden koşullarla birleştirdik. Çalışmada RSM ve ANFIS ile hesaplanmış iki sütun denetim hedefine alınmadı.
3. Kaynağın raporladığı logaritmayı `log10_eta_reported` alanında koruduk. `eta_mPa_s = 10**log10_eta_reported` ile fiziksel büyüklüğü geri elde ettik. Kaynak Şekil 8'de viskoziteyi mPa·s cinsinden gösterir. CSV'de dönüşümden sonra beş ondalığa yuvarlama, cihaz hassasiyeti belirtmez.
4. Yirmi deney 15 eşsiz kayma–derişim–tuz durumuna ve 9 eşsiz derişim–tuz formülasyonuna karşılık gelir. `group_id`, aynı formülasyondaki bütün kayma hızlarını ve tekrarları birlikte test dışında bırakmaya yarar. 30 s⁻¹, 2000 ppm, 12,5 g/L merkezinde altı tekrarın aritmetik ortalaması 7,00 mPa·s ve örnek standart sapması 0,32 mPa·s'dir. Yazarlar tekrarlar için ayrı numune lotları bildirmemiştir.
5. Şekil 8'de çalışma alanının dışındaki deneyler 5 kütle % NaCl ile gösterilir. Yoğunluk bilgisi olmadan bu değer eğitimdeki g/L birimiyle birebir eşleşmediği için bu deneyler sayısal tabloya eklenmedi. Çalışma içi gruplu sınama, başka bir laboratuvar veya polimer ailesindeki aktarım hatasını ölçmez.

Dosyalarda `partition=calibration`, gruplu çapraz doğrulamayı kurarken kullanılabilen kaynak satırlarını gösterir. Her dış döngüde bütün bir formülasyon kalibrasyonun dışında tutulur. Tablo 1'de her deney için tek kayma hızı seçilmiştir. Aynı ortalamaya sahip zincirler için genelleme ve sıfır-kayma viskozitesi ayrı deney gerektirir.
