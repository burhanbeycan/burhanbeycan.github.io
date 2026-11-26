# Seyreltik polimer çözeltilerinde viskozite

## Mol kütlesinden intrinsik viskoziteye

Model her çözücü ve sıcaklık için altı eşlenmiş polistiren fraksiyonuna `ln[η] = ln[η]ref + a ln(Mw/Mref)` bağıntısını eşit ağırlıkla uydurur. `[η]` mL/g, `Mw` g/mol cinsindendir. `Mref` eğitim kütlelerinin geometrik ortalamasıdır. Ampirik uyumun kütle ekseni, makalede bildirilen `Mw` değeridir. Geniş dağılımlı bir numunede viskozite ortalamalı `Mv` ile `Mw` farklı olabilir. Kayma hızı ve çözücü sıcaklığı ölçümler arası aktarıma etki edebilir. Veri kartındaki ± büyüklüklerinin istatistiksel düzeyi bilinmediğinden uyuma güven aralığı olarak girmezler. [Einaga vd., Tablo I](https://doi.org/10.1002/pol.1979.180171205)

| PS çözeltisi | Kalibrasyon, milyon g/mol | Üs | Fraksiyon tek tek ayrıldığında bağıl RMSE | Bağımsız ölçümde bağıl RMSE |
| --- | ---: | ---: | ---: | ---: |
| Benzen, 25,0 °C | 8,76–56,8 | 0,74764 | %2,554 (6) | %2,571 (1) |
| Siklohekzan, 34,5 °C | 8,76–56,8 | 0,50923 | %3,035 (6) | %16,265 (5) |

Tek-fraksiyon sınamasında eğitim ve sınama aynı makaledendir. JSON ayrıca her ana numuneden türeyen tüm fraksiyonları birlikte ayırır. Her çözücüdeki üç ek Einaga fraksiyonu ayrı bir düşük-kütle sınamasıdır. Bir yayını dışarıda bırakan sınamada NBS'nin benzen referansı 74,3 mL/g, model 72,39 mL/g verir. Baysal–Uyanık siklohekzan ölçümlerinin yalnız `10 000 000 g/mol` numunesi kalibrasyon aralığındadır. Bu numunede tahmin 281,92 mL/g, deney 211 mL/g'dır. İki ayrı sıcaklıktaki NBS ölçümü veri kaydında izlenebilir ve eş-sıcaklık skorunun dışında kalır. Buradaki bağıl RMSE, `100 √ortalama[((tahmin−ölçüm)/ölçüm)²]` değeridir. [NBS SRM 705](https://tsapps.nist.gov/srmext/certificates/archives/705_COA_1978-11-28.pdf) · [Baysal ve Uyanık, Tablolar 1–2](https://doi.org/10.1016/0032-3861(92)90695-S)

## Sonlu derişim ve kapiler sınama

`C = c[g/L]/1000` g/mL ve `u = C[η]` boyutsuzdur. Huggins açılımı `ηr = 1 + u + kH u²`, `ηsp = ηr−1` verir. Aynı akış durumundaki saf çözücünün dinamik viskozitesi biliniyorsa `ηçözelti = ηçözücü ηr` elde edilir. `cη = 1000/[η]` g/L viskometrik bir derişim ölçeğidir. Önceki uygulamada `Rg` ile bulunan geometrik örtüşme derişimi için ayrıca zincir konformasyonu ve hacim tanımı gerekir. `u ≪ 1` seyreltik açılımın doğal sınırıdır. `u = 1`, tüm numuneler için deneysel bir keskin geçiş derişimi sağlamaz.

Huggins için `ηsp/C` ve Kraemer için `ln(ηr)/C` ayrı ayrı derişime doğrusal uydurulur. Sıfır-derişim kesişimleri ayrı raporlanır. Kraemer katsayısı işaretli biçimde `ln(ηr)/C = [η]K + kK[η]K² C` tanımlanır. Aynı sıfır-derişim limitine sahip iki asimptotik açılım `kK = kH − ½` verir. Sonlu aralıkta ayrı regresyonların bu eşitliği tam sağlaması beklenmez. Kaynakta ilgili fraksiyon için `kH` eksikse sıfır derişim limiti ile viskometrik ölçek gösterilebilir.

| Seri | Sınama c, g/L | Eğitimden bulunan u | Ölçülen ηr | Huggins ηr |
| --- | ---: | ---: | ---: | ---: |
| PS/toluen, Copenhaver vd., 30 ± 1 °C | 19,930 | 1,360 | 3,6300 | 3,3612 |
| PS/toluen, PHYWE, 25 °C | 20,000 | 1,302 | 2,8890 | 2,8446 |
| PS/siklohekzan, Luh, 633 000 g/mol | 9,956 | 0,662 | 1,9161 | 1,9126 |
| PS/siklohekzan, Luh, 1 770 000 g/mol | 4,704 | 0,521 | 1,6821 | 1,6826 |

Her satırda aynı numunenin en yüksek derişimi uyumdan çıkarıldı. Toluendeki iki yüksek-derişim noktası seyreltik açılımı zorlar. Copenhaver serisinde farklı öğrenci ekipleri farklı viskozimetreler kullandı. Kaynakta sorgulanan ikinci gözlem hesapta tutuldu. Her iki toluen serisinde ölçülen sürelerden kaynağın yayımladığı boyutsuz bağıl viskozite kullanıldı. Luh çalışmasında yoğunluk ve kalibrasyon kullanılarak hesaplanmış dinamik viskoziteler de tabloda mevcuttur. Dört seri farklı cihaz, numune ve sıcaklıklardadır. Kaynakta aynı beş derişimden türetilmiş `[η]` değerleri, tüm seriyle uyumun yeniden üretilebilirliğini sınar. [Copenhaver vd., Tablo 1](https://doi.org/10.1021/acs.jchemed.4c00439) · [PHYWE, Tablo 1](https://spegroup.ru/upload/wikifiles%202/Experiment_guide_P3010601.pdf) · [Hwa Luh, Tablo 8](https://scholarsmine.mst.edu/masters_theses/6870/)
