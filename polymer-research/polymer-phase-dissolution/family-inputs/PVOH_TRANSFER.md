# PVOH’da çözücüye özgü öğrenme

## Soru ve hedef

Aynı PVOH numunesinin su ve DMSO’daki `[η]` değerlerinin her ikisini de sakladığımızda, mol kütlesinin ve hidroliz derecesinin etkilerini çözücüler arasında kısmen paylaşmak hatayı azaltıyor mu? Hedef seyreltik sınırdaki intrinsik viskozitedir (mL/g). Hata ölçümleri her koşulda önce mL/g uzayında gözlem ve tahmin farkından hesaplanır. MAPE her satırın mutlak bağıl hatasının ortalamasıdır.

## Öğrenilen modeller

Her çözücü için yerel model `ln([η]/(1 mL/g)) = b_s + a_s ln(Mw/(50 000 g/mol)) + h_s (DH − 80)/10` biçimindedir. DH molce yüzde puan cinsindendir. Pozitif sonuç logaritmik yanıtla sağlanır. İntrinsik viskozite, Huggins katsayısı, örtüşme derişimi, hedef çalışmanın ölçülen aynı numunedeki diğer çözücü sonucu veya [η]'den türetilmiş mol kütlesi girdi olarak kullanılmaz.

1. **Mol kütlesi.** Önceki sürümün ortak mol kütlesi üssü ve çözücü başına ayrı kesişimi.
2. **Ortak hidroliz etkisi.** Önceki sürümün her iki çözücü için ortak kütle ve hidroliz eğimleri. Kesişimler yine çözücüye özgüdür.
3. **Birincil aday.** Ortak kütle üssü korunur. Su–DMSO arasındaki eşlenmiş log-viskozite farkının yarısı, eğitim numunelerinin DH değerine karşı uyarlanır. Eğimin yarısı pozitif/negatif işaretle ortak hidroliz eğimine eklenir. Su ve DMSO sonuçlarının aynı numuneye ait olmasından faydalanır. Bu birincil aday sonuçlar hesaplanmadan belirlenmiştir. Önceki beş numunedeki temel model sonuçları önceden biliniyordu.
4. **İki eğimde kısmi havuzlama.** Daha sonra, kütle üssü ve DH eğiminin ikisi için çözücüye özgü katsayılar öğrenilmiştir. Her katsayı yerel ve ortak katsayının aritmetik ortalamasıdır. Kesişim her çözücü için merkezlenmiş eğitim verisinden hesaplanır. Bu aday keşifseldir. Aynı beş numunenin ilk model değerlendirmesi görüldükten sonra karşılaştırmaya eklenmiştir.
5. **Sulu kaynaktan aktarım.** Nagy'nin 12 sulu SEC ölçümünde iki kaynak eğimi öğrenilir. Suda bu iki eğim dördüncü modeldeki su eğimleriyle eşit ağırlıkla birleştirilir. Böylece kaynak ortamının kendine özgü kesişimi aktarılmaz. DMSO kolunda model 4'ün tahminleri geçerlidir. Aktarımın dışarıdan veri sağlaması, sıcaklık ve ortam farklarını ortadan kaldırmaz.

Yarıya küçültme sabittir, dış hedef hatasına göre ayarlanmaz. Katsayılar yalnız eğitimdeki dereceler üzerinden hesaplanır. Her saklı numune için iki gözlem birlikte dışarıda kalır. Dört hedef-eğitim numunesinde sekiz yanıt vardır. Yeni tahmin için beş Parisi numunesiyle tam uyum yapılır. Beş numune üzerinden kalibre edilmiş bir tahmin aralığı sunulmaz.

## Eğitim ve bağımsız seriler

| Kaynak | Gözlem | Sınamadaki görev |
|---|---:|---|
| Parisi ve ark., *Journal of Rheology* 66, 1141–1150 (2022), Tablolar I ve III | 5 PVOH × su/DMSO = 10 | 25 °C hedef serisi. Boyut ve kimya önceki PVOH veri notuyla aynıdır. |
| D. J. Nagy, *Handbook of Size Exclusion Chromatography and Related Techniques*, bölüm 10 (2004), Tablolar 1 ve 3. 2003 *American Laboratory* çalışmasının aynı tablosu | 12 | Yalnız açık sayısal hidrolizli sınıflar. %99+ için tam değer belirtilmediği için dışarıda. Kalın sütun MALLS M𝑤, [η] farklı viskozimetre sinyalidir. Tablo 3 dL/g değerleri 100 ile çarpıldı. 35 °C, 0,05 M NaNO₃ sulu elüent. Kaynak yalnız beşinci modelin su eğimlerini etkiler. |
| Gebben, van den Berg, Bargeman ve Smolders, *Polymer* 26, 1737–1740 (1985), Tablo 1, Tablo 2, Şekil 1 | 2 | Çapraz bağlama öncesi başlangıç numuneleri. M𝑤 HPLC–düşük açılı ışık saçılması, DH sabunlaşma, [η] kılcal viskozimetri. 25 °C, PVA-126-98 için Tablo 2 ile doğrulanıyor. PVA-96-88’de Tablo 1’de sıcaklık ayrı verilmemiştir, Şekil 1 yaklaşık 25 °C başlangıç noktasına işaret eder. Çapraz bağlanmış fraksiyonlar ve 40 °C değerleri dışarıda. Her iki DH değeri hedef eğitim desteğinin üzerindedir. |
| Korelc ve ark., *European Journal of Pharmaceutical Sciences* 192, 106645 (2024), Tablo 2, bölüm 2.2.1.3 | 2 | MALLS ve viskozimetre dedektörlü SEC. Ortam 10 mM NaCl ve %0,002 NaN₃ içeren su. SEC sıcaklığı belirtilmez. DLS veya Brookfield deneylerindeki 25 °C SEC için kullanılamaz. PVA35, PVA50 ve PVA74 numunelerinde [η] belirlenememiştir. Yalnız PVA88 ve PVA98 aktarım kontrolüne girer. |

Gebben kapsamındaki en kesin sıcaklık eşleşmesi, PVA-126-98 numunesinde 105,1 mL/g'dır. Ortak etki modeli 88,30, iki eğimli model 96,71 ve sulu kaynak aktarımlı model 81,45 mL/g tahmin eder. Böylece yardımcı kaynaktan aktarımın bağımsız bir ölçümde hatayı büyütebildiği görülür. İki farklı laboratuvarın sonuçları tek bir RMSE altında toplanmamıştır.

## Doğrulama sonuçları

| Sınama | Mol kütlesi | Ortak DH | Çözücü DH | İki eğim | Sulu aktarım |
|---|---:|---:|---:|---:|---:|
| Parisi, su ve DMSO, 10 koşul | 13,64 | 13,51 | 12,75 | 10,21 | 10,11 |
| Parisi, su, 5 | 13,87 | **6,34** | 6,45 | 7,67 | 7,40 |
| Parisi, DMSO, 5 | 13,41 | 18,03 | 16,84 | **12,24** | **12,24** |
| Gebben, 2 başlangıç numunesi | 37,80 | 17,68 | 13,13 | **12,82** | 22,23 |
| Korelc, 2 farklı koşullu ölçüm | 19,63 | 11,83 | 10,48 | **9,19** | 11,20 |

Tabloda RMSE birimi mL/g'dır. Suda en iyi önceki model ortak hidroliz katsayısıyla 6,34, DMSO'da ise yalnız mol kütlesiyle 13,41 mL/g verir. Aynı çözücüdeki koşullar arasında daha iyi sonucun seçilmesi ve aynı veride puanlanması yeni ve doğrulanmış bir seçim kuralı oluşturmaz. Birleşik puan çözücü başına beş ölçümden hesaplanmıştır. DMSO’da mutlak viskozite daha yüksek olduğu için aynı bağıl sapma birleşik RMSE’ye daha fazla katkı verebilir. MAPE ve MAE arayüzde aynı satırlardan hesaplanır.

[Parisi 2022](https://doi.org/10.1122/8.0000435) · [Nagy, s. 272 ve 277](https://www.marinalg.org/wp-content/uploads/2012/09/Wu-2004-Handbook-of-size-exclusion-ch10.pdf) · [Gebben, tablolar](https://ris.utwente.nl/ws/portalfiles/portal/256201525/Gebben1985intramolecular.pdf) · [Korelc, Tablo 2 ve SEC yöntemi](https://research.chalmers.se/publication/538629/file/538629_Fulltext.pdf)

Hesapları yinelemek için proje kökünde `node scripts/build-pvoh-transfer.cjs` ve `node tests/test-pvoh-transfer.cjs` kullanılır. Veri ve sonuç dosyaları çevrimdışı da kullanılabilir.
