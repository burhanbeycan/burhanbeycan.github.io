# Ölçülen zincir boyutundan intrinsik viskoziteye aileler arası aktarım

## Eşlenmiş deneyler

Fetters, Hadjichristidis, Lindner ve Mays'ın kritik derlemesinin 2, 4 ve 6. tablolarından aynı satırında M, [η] ve Rg birlikte bildirilen lineer homopolimerleri aldık. Bunlar eğrilerden tahmin edilmiş sanal noktalar değildir. Siklohekzan iyi çözücüdür. Derlemedeki koşul 25 °C'dir. Yazarlar yakın sıcaklıklardaki iyi çözücü ölçümlerini birlikte tabloya aldıklarını bildirir. CSV’de 298,15 K bu tablolardaki nominal koşuldur. Her kaynağın termostat sıcaklığının satır bazında bağımsız doğrulandığı anlamına gelmez.

| Polimer | Seçilen satırlar | Birincil deney kaynağı |
|---|---:|---|
| 1,4-Polibütadien | 4 | Tablo 2, kaynak 24 · [Roovers ve Martin](https://doi.org/10.1002/polb.1989.090271209). Makaledeki lineer numuneleri içeren derleme satırları kullanıldı. |
| 1,4-Poliizopren | 6 | Tablo 4, kaynak 30 · [Tsunashima ve ark.](https://doi.org/10.1021/ma00182a043), kaynak 31 · [Davidson ve ark.](https://doi.org/10.1021/ma00176a048). |
| Poliizobütilen | 12 | Tablo 6, kaynak 32 · [Matsumoto ve ark.](https://doi.org/10.1002/pol.1972.160100103), kaynak 36 · [Nakamura ve ark.](https://doi.org/10.1021/ma00029a020). [Nakamura'nın deneylerinin açıklaması](https://ir.library.osaka-u.ac.jp/repo/ouka/all/27674/10137_Dissertation.pdf). |

Tablo 6’daki kaynak 34’ün iki tam satırı çıkarıldı. Birincil yayında numune *butil kauçuk* olarak tanımlanır. [Tsuji ve Fujita](https://doi.org/10.1295/polymj.4.409). PB için 1. tablonun 25 numaralı kaynağında bazı θ yarıçapları makaledeki farklı büyüklükten hesaplanmıştır. Bu tablo ve koşul bildirilen siklohekzan verisine dahil edilmedi. Derlemenin M sütunundaki sayısal değerlere sadık kalındı. Ortalaması birincil kaynakta belirlenebilen poliizobütilen serileri Mw olarak, kalan seriler M_table olarak işaretlendi. Tablolar dar dağılımlı fraksiyonları içerir. [Fetters ve ark., tablolara açık erişim](https://srd.nist.gov/jpcrdreprint/1.555949.pdf).

## Model ve birimler

Ölçülen statik zincir boyutu Rg ile viskozimetrinin etkin hacmini Flory–Fox ölçeklenmesiyle ilişkilendiriyoruz. Bir molekül için Vg = 4π Rg³/3 eşdeğer geometrik hacimdir. Boyutsuz U = [η] M /(Nₐ Vg) tanımında M g/mol, Rg nm, [η] mL/g olarak alınır; 1 nm³ = 10⁻²¹ mL ve Nₐ = 6,02214076 × 10²³ mol⁻¹. Tahmin [η] = U × Nₐ Vg/M olur. U her bölümde yalnız eğitim numunelerinden öğrenilir. Bu fiziksel katsayı, çözücü kalitesi, zincir konformasyonu ve kullanılan dağılım ortalamalarının etkisini taşıyabilir. [Viskozite yarıçapı ve evrensel oranın tanımı](https://doi.org/10.1021/ma500960f).

Kütle modeli ln[η] = b + a ln M şeklinde uydurulur. İki model tam olarak aynı numunelerle kurulur ve aynı saklı hedefleri tahmin eder. U için logaritmik ortalama, kütle için logaritmik doğrusal regresyon kullanılır. Her eğitim kimyası eşit toplam ağırlık alır ve bu ağırlık o kimyanın fraksiyonları arasında eşit bölünür. Boyut modelinin ayrıca **hedef numunenin deneysel Rg değerini** kullandığı özellikle önemlidir. Bir giriş yapısından Rg üreten model bu sınamada değerlendirilmedi.

## Bölmeler ve yorum

- **Tek kimya.** Üç polimer sırayla tümüyle ayrılır. İki polidienin kimyasal yakınlığı korunur.
- **Omurga.** PB ve PI birlikte dışarıda bırakılır. Ters yönde iki polidienle kurulan model doymuş omurgalı PIB'nin tüm kayıtlarını tahmin eder. Bu ters yön, PIB'nin tek kimya bölmesiyle aynıdır.
- **Harici kimya ve ortam.** 22 siklohekzan gözleminden elde edilen iki model aynen korunur. Birinci sınamada, 25 °C’de benzendeki altı PS fraksiyonunun [η] ve Rg değerleri önceki Studio modülündeki fraksiyon kimliğiyle eşleştirilir. İkinci sınamada 25 °C’de suda kütlece %0,2 sodyum azit içeren ortamda ölçülen yedi nötr, lineer pullulan fraksiyonu kullanılır. Orijinal tabloda şüpheli görülen iki ek satır ana sınamaya girmez. Her iki sınamada aynı anda çözücü ve kimya değişir. [PS boyutları](https://doi.org/10.1021/ma60066a022) · [PS viskozimetrisi](https://doi.org/10.1002/pol.1979.180171205) · [Pullulan birincil verileri](https://doi.org/10.1021/ma00139a016).
- **Hata.** Her bölümde MAPE = 100 × ortalama(|tahmin − ölçüm|/ölçüm). MAE mL/g cinsindendir. Bağıl RMSE ve işaretli ortalama yüzde hata da verilir. Tüm kayıtlar saklı hedef için aynı ağırlıktadır. Aralık etiketleri yalnız eğitimde görülen en küçük ve en büyük sayısal değerleri karşılaştırır.

Boyut modelinin girdisi Rg’nin küpünü içerir. Örneğin M ve U sabitken Rg'nin %5 büyük ölçülmesi tahmini (1,05³ − 1) × 100 = %15,76 büyütür. Bu hesap bir duyarlılık örneğidir. Kaynaklar her fraksiyon için deneysel standart hata vermediğinden model tahminlerine bu oranla güven aralığı atanmamıştır. Sıcaklık ve çözücü kalitesi sabit kabul edildiğinde bile viskozimetri ve ışık saçılmasının polidispersiteye duyarlılığı farklıdır.

## Yeniden üretim

Site düzeninde `node scripts/build-hydro-family.cjs dist` ve `node tests/test-hydro-family.cjs` komutlarını kullanın. Düz GitHub Pages veya çevrimdışı paket için ilk komutta `dist` yerine `.` verin. Tüm kaynak satırları, saklı numune kimlikleri, katsayılar ve satır bazında tahminler çıktı dosyasında bulunur.
