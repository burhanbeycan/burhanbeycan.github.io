# Huggins denklemiyle sınırlandırılmış Gauss süreci

Polimer viskozitesi literatüründeki fizik-kısıtlı öğrenme yaklaşımı, bilinen bağıntıyı modelin içine yerleştirir. PENN çalışması saf **polimer eriyiklerinde** sıcaklık, mol kütlesi ve kayma hızının etkilerini inceler. Buradaki deney ise **seyreltik PS çözeltilerinin** derişime bağlı bağıl viskozitesini, kapiler verilerle sınar. Eriyik çalışmasından taşınan öğe fiziksel yasayla koşullandırma fikridir. Model ve parametreler bu çözelti problemine özgüdür. [PENN çalışması](https://doi.org/10.1038/s41524-025-01532-6)

## Model ve kalibrasyon

`u = c[g/L] [η][mL/g]/1000` alındığında `ηr = 1 + u + kH u² + δ(u)` kullanılır. Artık terim sıfır ortalamalı Gauss süreci ile öğrenilir. Kovaryans `K(u,v) = σ² u² v² exp[−(u−v)²/(2ℓ²)]` olur. Çarpanlar hem `δ(0)=0` hem `δ′(0)=0` koşulunu sağlar. Saf çözücü limitinde `ηr=1` kalır, başlangıç eğimi `[η]` korunur. Referans koşul `ℓ=0,5`, `σ=0,25` ve ηr biriminde düzenlileştirme standart sapması `0,02` kullanır. Bu üç sayı varsayımdır. Posterior çözüm Cholesky ayrıştırmasıyla hesaplanır. Koşullu standart sapma sabit Huggins ortalaması ve çekirdek altında gizli artık fonksiyona aittir. Deney hatası, `[η]` önkestirimi, seçilen fonksiyon ve çekirdek parametrelerindeki belirsizlikler ayrıca değerlendirilmelidir. [Rasmussen ve Williams, bölüm 2](https://gaussianprocess.org/gpml/chapters/RW2.pdf)

**Numune içinde** Huggins ve GP aynı eğitim noktalarını görür. `[η]` ve `kH` yalnız bu noktalarda Huggins regresyonuyla hesaplanır. Üç tekrar kullanılabilir bölme vardır. En yüksek derişim, en yüksek iki derişim veya üç en küçük derişimden ortadaki ayrılır. Sonuncusunda iki eğitim noktası Huggins doğrusunu tam belirlediği için öğrenilecek artık kalmaz ve iki model çakışır. Farklı bölmelerin sınama gözlemleri kesişir, toplamları bağımsız sınama sayısı vermez.

**Bir PS fraksiyonu tamamen dışarıdayken** iki molekül ağırlığındaki `[η]` için diğer yayındaki altı yüksek-kütle fraksiyonundan, aynı çözücü/sıcaklıkta Mark–Houwink önkestirimi kullanılır. Eğitim fraksiyonunda bu önkestirim sabit tutulur. `kH = Σu²(ηr−1−u)/Σu⁴` ile bağıl viskozitede kare hatası en aza indirilir. Artık aynı eğitim noktalarında öğrenilir. Bu `kH`, `[η]` önkestirimindeki sapmayı da yansıtabilecek bir katsayıdır. Sınama fraksiyonunun viskozite tablosu veya bu tablodan türetilen `[η]` eğitime girmez. Eğitim ve sınama fraksiyonlarının Luh'un aynı çalışmasından gelmesi numuneler arası sınama sağlar. Her iki mol kütlesi, bağımsız `[η]` önkestiriminin `8,76–56,8 milyon g/mol` aralığının altındadır. [Luh, Tablo 8](https://scholarsmine.mst.edu/masters_theses/6870/) · [Einaga vd., Tablo I](https://doi.org/10.1002/pol.1979.180171205)

| Sınama | Ayrılan ölçüm | Huggins bağıl RMSE | Referans GP bağıl RMSE | 9 ayarda GP bağıl RMSE aralığı |
| --- | ---: | ---: | ---: | ---: |
| 633 000 g/mol, tüm fraksiyon | 5 | %0,604 | %0,411 | %0,277–0,582 |
| 1 770 000 g/mol, tüm fraksiyon | 5 | %0,585 | %0,382 | %0,327–0,541 |
| Copenhaver PS/toluen, en yüksek c | 1 | %7,406 | %10,244 | %0,982–21,974 |
| PHYWE PS/toluen, en yüksek c | 1 | %1,537 | %1,943 | %1,681–2,120 |

Tüm-fraksiyon sınamasında iki hedef için de iyileşme görüldü. İki toluen serisinin en yüksek derişiminde referans GP hata artırdı. Bu iki seride `u>1` ve sınama noktası eğitim derişimlerinin üstündedir. Dokuzlu duyarlılık incelemesi, referans sonuçlar görüldükten sonra keşifsel olarak belirlendi. `ℓ ∈ {0,25; 0,5; 1}` ile `σ ∈ {0,125; 0,25; 0,5}` kombinasyonları, sabit `0,02` düzenlileştirme ile tümüyle raporlanır. Tüm karşılaştırmalarda eğitim ve sınama kümeleri aynıdır. Altı sınamanın her ayarı JSON'da bulunur. Referans hiperparametreler değiştirilmez.

Arayüzdeki değerlendirme alanı, derişim (`u`) aralığını ve `[η]` için mol kütlesi kalibrasyonunu ayrı ayrı izler. `u<1`, gözlenen eğitim aralığı ve pozitif/artan yanıt bir ön kontroldür. Koşullu bant için deneysel kapsama oranı kalibre edilmemiştir.
