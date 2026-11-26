# PVA–su için akış parametrelerinin ters kestirimi

Carrasco-Venegas ve arkadaşları, Peru'da dağıtılan Merck 141351 PVA ürünü için kütlece %4, %6, %8 ve %10 derişimde, 20, 25, 30 ve 35 °C'de reometri ölçümleri yayımladı. Bu klasörün hedefi bu ürüne ait kararlı kayma viskozitesidir. Numunenin sayısal Mₙ, M𝑤, Đ ve hidroliz derecesi makalede bildirilmemiştir. [Birincil makale, §2.3 ve Tablo 2](https://doi.org/10.3390/polym15071743).

## Veriden fiziksel modele

1. `pva_table2.csv`, birincil makalenin Tablo 2'sindeki 244 pozitif deneysel gerilim–hız çiftini içerir. Konsantrasyon kütlece yüzdedir. Her gözlem için η [cP] = 1000 × τ [Pa] / γ̇ [s⁻¹] hesaplanır. τ yalnızca gözlem etiketini oluşturmakta kullanılır.
2. Her derişim–sıcaklık eğrisinde, x = ln(γ̇ / 30 s⁻¹) olmak üzere ln(η/cP) = A − b₁ min(x,0) − b₂ max(x,0) bağıntısı tersine çözülür. A, 30 s⁻¹'deki viskozitenin logaritmasıdır. 0 ≤ bᵢ ≤ 1 kısıtları pozitif viskoziteyi, kayma incelmesini ve azalmayan kayma gerilimini korur. b₁ ve b₂ farkına verilen küçük karesel ceza seyrek örneklenen eğrileri kararlılaştırır. 30 s⁻¹ sabit sayısal düğümdür. Bu değere malzemeye özgü bir gevşeme frekansı anlamı verilmez. Tasarım matrislerinin koşul sayısı ve her eğri için hata `inverse_parameters.csv` içindedir.
3. `inverse_parametric` yöntemi, 16 eğrinin her birine eşit ağırlık vererek üç ters kestirim parametresinin derişim/sıcaklık bağımlılığını Ridge regresyonu ile öğrenir. Girdi tabanı [1, u, v, u², uv] olup u=ln(c / 7 wt%), v=1000(1/T − 1/298,15)/0,08. Kestirimde eğimler fiziksel aralığa sınırlandırılır. Tahmin işlevi yalnız eğitim derişim/sıcaklık aralığında ve komşu hücrelerin hepsinin ölçülen kayma penceresinde değer verir. η₀, sonsuz kayma platosu, aktivasyon enerjisi veya kinetik gevşeme süresi bu katsayılardan çıkarılamaz.
4. `shared_power_law` bütün koşullarda tek kayma üssüyle basit fiziksel karşılaştırmadır. ln(c) ve 1/T kaymaları içerir. `direct_ridge` aynı veriye uydurulan, x² ve etkileşimleri de içeren serbest log-viskozite regresyonudur. `inverse_curve_ml`, tersine bulunan üç parametrenin ortak kuvvet yasasından sapmasını öğrenen Gauss çekirdekli Ridge modelidir. Ana model diferansiyel denklem veya moleküler-dinamik vekili içermez.
5. Kullanıcı aynı formülasyona ait tek bir referans η ölçümüne sahipse yalnız genlik bu ölçümün tahmine oranıyla yeniden ölçeklenebilir. Çalışma sırasında yeni bir ölçüm istenmez. Bu koşullu sonuç deneysel desteği olmayan yeni bir numunede başarıyı göstermez.

## Dahili hata ve bağımsız ürün sınaması

Her dış kat için tüm parametreler yalnız eğitim kümesinde yeniden kestirildi. Aşağıdaki MAPE (%) her bağımsız akış eğrisinin ortalama mutlak yüzde hatasını eşit ağırlıkla birleştirir. Her protokolde 16 eğrinin bütün 244 hedefi bir kez testte yer aldı. Uç derişim veya sıcaklık dışarıdayken yapılan ekstrapolasyonlar zorluğu görünür kılmak için benchmark'a dahil edildi. Uygulama komşu ölçüm penceresi dışına çıkan bir kullanıcı isteğine sayısal sonuç vermez.

| Model | Tüm derişim dışarıda | Tüm sıcaklık dışarıda | Tüm derişim–sıcaklık eğrisi dışarıda |
|---|---:|---:|---:|
| Ortak fiziksel kuvvet yasası | 38,32 | 19,85 | 23,27 |
| Doğrudan Ridge | 19,05 | 7,13 | 8,37 |
| Ters eğri + çekirdek düzeltmesi | 27,32 | 8,70 | 11,47 |
| **Ters eğri + koşula bağlı parametre regresyonu** | **18,13** | **6,04** | **6,67** |

Tek referans sınamasında her ayrılmış eğrinin 30 s⁻¹'e en yakın ölçümü yalnız kalibrasyona ayrıldı. 16 referans çıkarıldıktan sonra aynı 228 noktada tahminler karşılaştırıldı. Eğri dışarıda bırakıldığında eşit eğri ağırlıklı MAPE %6,66'dan %3,00'a, tüm derişim dışarıda bırakıldığında %18,09'dan %5,28'e değişti. Aynı laboratuvar ve aynı ürüne ait bu yöntem geliştirme sonuçları bağımsız çalışma doğrulaması sağlamaz. `benchmark.json` ve `one_reference_benchmark.json` ham nokta ağırlıklı metrikleri de içerir. `equal_curve_mean.median_ape_pct`, eğri içi medyanların ortalamasıdır.

Dış sınamada dondurulmuş Carrasco modeli Pérez-Puyana ve arkadaşlarının **başka bir PVA ürününe**, 298,15 K ve 10 s⁻¹'de uygulandı. Yazarlar M𝑤 = 130.000 g/mol ve %86,7–88,7 hidroliz bildiriyor. Derişimler kütle esaslı ve çözücü sudur. %5, %7,5 ve %10 için gözlenen viskoziteler sırasıyla 126, 405 ve 1490 cP, tahminler 131,68, 170,29 ve 249,03 cP'dir. Üç noktada MAPE %48,58, ±%10 içindeki nokta sayısı 1'dir. Ürün, laboratuvar ve çözelti hazırlama prosedürü birlikte değiştiği için sapmanın tek bir etkene atfedilmesi mümkün değildir. %2,5 ölçümü ve Turan ve arkadaşlarının beş noktası sayısal destek dışında kaldı. [Pérez-Puyana, Tablo 1](https://doi.org/10.1088/2053-1591/aab164) · [Turan, Tablo 3](https://doi.org/10.35378/gujs.794638).

## Çalıştırma

Python 3, NumPy ve SciPy gerekir. Tarayıcıdaki hesap için yalnız `inverse-engine.js` ve `model_state.json` gerekir. Dosyalardaki nihai modeller bütün geliştirme verileriyle eğitilmiştir. Tam hata tabloları dış tahminlerden hesaplanır.

```bash
python pva_model.py build
python pva_model.py predict --wt 4 --kelvin 298.00 --rate 10
python pva_model.py predict --wt 7.5 --kelvin 298.00 --rate 10 --anchor-rate 30 --anchor-cp 80
python external_product_test.py --measurements ../literature/pva_measurements.csv
node test-parity.cjs
```

`build` hedef ve kat düzeyindeki çıktıları yeniden oluşturur. `heldout_predictions.csv`, `one_reference_predictions.csv` ve `fold_metrics.csv` her dış kat için raporlama yapılmasını sağlar. `external_product_test.py` açıkça tanımlanan konsantrasyon, sıcaklık, kayma hızı, hedef ve çözücü ölçütlerine uyan dış satırları seçer ve dışarıda kalan koşulları kaydeder. JS arayüzünde `PvaInverseModel.predict(state, wtPercent, temperatureK, shearRate)`, `eta_cP`, `supported` ve `reason` alanlarını döndürür. Tarayıcı ile Python arasında 288 noktada sayısal eşdeğerlik ve koşul desteği kontrolü vardır.

Ek dış ürün deneyi aynı Merck modelini Parisi vd. Şekil 4'ün sulu çözelti ölçümlerine uyguladı. Kaynak şeklin sayısallaştırmasında dört kimliği doğrulanmış ürün kaydedildi. 283 noktanın 35'i sayısal çalışma pencereleriyle uyuştu (kütlece %4, 298,15 K). Bu 35 noktanın hiçbiri ±%10 içinde kalmadı. Ürün başına MAPE 144k–88 için %245,1 (10 nokta), 171k–80 için %144,3 (10), 24k–80 için %1.988,2 (5), 49k–75 için %1.714,8 (10) oldu. Küçük ölçülen viskoziteler yüzde hatayı büyüttüğünden dosyada logaritmik ve cP ölçeğinde hatalar da bulunur. Bu karşılaştırmada Parisi verisinden hiçbir katsayı öğrenilmedi. [Şekil 4 ve deney yöntemi](https://doi.org/10.1122/8.0000435).
