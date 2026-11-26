# Faz Dengesi ve Çözünme

Polymer Research Studio’nun polimer–çözücü uygulaması. Tarayıcıda hesaplar. Ana sayfada faz diyagramı, ara yüz dinamiği, atomistik film ve yapı aktarımı açılır.

## Çalıştırma

`index.html` ile aynı dizinde `python -m http.server 8000` çalıştırıp tarayıcıda `http://localhost:8000` adresini açın. Sites sürümünde web dosyaları `dist/` içindedir. GitHub Pages inceleme sürümünde doğrudan proje dizinindedir.

## Model ve gözlenebilir

| Ekran | Gözlenebilir ve kapsam |
| --- | --- |
| PS–siklohekzan | Bulutlanma sıcaklığına uyarlanan sıkıştırılamaz, ikili Flory–Huggins modeli. İki model karşılaştırılır. χ = A + B/T ve χ = A + B/T + Cφ. Bileşim polimer hacim kesridir. 0,10 MPa. Dar dağılımlar için tek etkin zincir uzunluğu kullanılır. Denge fazları ortak teğetle, spinodal türevle belirlenir. |
| PEtOx–su | Numuneye özgü optik bulutlanma ölçümleri ve sınırlı derişim aralığında interpolasyon. Kalorimetrik geçiş, optik bulutlanma ve moleküler hidratasyon ayrı gözlenebilirlerdir. |
| İki boyutlu alan | Korunumlu Cahn–Hilliard dinamiği, periyodik sınırlar, sabit hareketlilik ve gradyan katsayısı. Zaman ile mesafe indirgenmiştir. Homojen sıçrama ve çözücüyle temas eden tabaka ayrı başlangıç koşullarıdır. |
| Atomistik film | 298,15 K ve 1 bar’da 700 TIP3P su molekülü içindeki iki H/H uçlu PEtOx trimeri. 8 ps kısıtlı hazırlığın ardından iki ayrı 60 ps yörünge. Film, ilk yerel su girişi etrafındaki kareleri gösterir. Etkileşimli oynatıcıda her iki yörüngenin tamamı vardır. |
| Yapı aktarımı | Önceki modülün yerel JSON çıktısını veya kimyasal editörü kullanır. Tek bir oligomerin formülü, stereokimyası, uçları, koordinatları ve mol kütlesi bir deneysel numunenin mol kütlesi dağılımından ayrı tutulur. |

## Denetim ve sonuç

- N = 100, χ = 0,65 ve φ = 0,020 kontrolü pozitif eğrilik ve metastabil homojen durum verir. Denge bileşimleri 0,009768 ve 0,271968, polimerce zengin fazın hacim payı %3,90’dır.
- PS verisi 3PDB’de birincil yayın kimliğiyle süzülür. Siporska ve arkadaşlarının 73 değeriyle uyarlanan bileşime bağlı model, Saeki ve arkadaşlarının ayrı 24 bulutlanma noktasında 0,77 K RMSE verir. Aynı sınama kümesinde bileşimden bağımsız model 2,45 K RMSE verir. Üç yeni mol kütlesini ayrı seçmek mümkündür.
- Özgün yayın, numune, kompozisyon temeli ve model tahmini her bulutlanma kaydıyla birlikte görüntülenebilir. Eşfaz bileşimlerinin ve spinodalin bağımsız deneysel sınaması açık araştırma sorusudur.
- PEtOx eğrisinde Aquazol 50'nin deneysel derişim birimi g/L olarak korunur. Farklı laboratuvarların ısıtma hızları ve geçiş tanımları gösterilir. Kalorimetrik büküm noktaları ayrı tabloda kalır.
- Kayıtlı yörüngeler OpenMM GAFF2/AM1-BCC ve TIP3P hesaplarından gelir. Filmin oynatma süresi ile MD saati ayrıdır. Yerel su girişi ve hidrojen bağı geometrisi bu iki kısa koşudan okunabilir. Filmden kütlesel çözünme süresi türetilemez.

## Kaynaklar

- [3PDB bulutlanma kayıtları](https://pppdb.uchicago.edu/cloud_points)
- [Siporska ve ark.](https://doi.org/10.1039/B301190A) ve [Saeki ve ark.](https://doi.org/10.1021/ma60032a020)
- [Hacioglu ve ark., PEtOx ölçümleri](https://doi.org/10.1021/acs.macromol.4c01341) ve [Zhang, Tablo 2-2](https://backoffice.biblio.ugent.be/download/5878271/5878291)
- [Gauche ve Felisberti, kalorimetrik geçişler](https://doi.org/10.1021/acsomega.9b01269)
- [Orijinal atomistik yörüngeler ve protokoller](https://github.com/burhanbeycan/burhanbeycan.github.io/tree/ca4c599eb76017723de36b225b0c793d32e7deb6/polymer-research/polymer-solvent-md/hydration/data/pair)

## Yeniden üretim

`node scripts/test_phase.cjs`, `node tests/kinetics.test.cjs` ve `python tests/check_distribution.py` dağıtım ile sayısal sonuçları sınar. PS uyumunu yeniden kuran `scripts/fit_ps.py` için Python, NumPy ve SciPy gerekir. `python scripts/render_contact.py --output ./film` komşu `polymer-solvent-md/hydration/data/pair` dizinindeki yörüngeleri okur. Özel konum ve seçenekler `--help` ile listelenir. Çizim için NumPy, Pillow, FFmpeg (libx264) ve DejaVu Sans gerekir. OpenChemLib JS sürümü ve BSD lisansı dağıtılan dosyanın başında bulunur.

## Ozmotik Basınç ve Çözelti Etkileşimleri

Bağımsız giriş `osmotic.html` üzerinden açılır. PS–siklohekzan faz ekranı numuneyi, integral χ modelini, sıcaklığı ve polimer hacim kesrini bu uygulamaya aktarır. Kütlesel derişim çözeltinin litresi başına gramdır. Saf çözücünün referans molar hacmi 108,04 mL/mol, PS özgül hacmi 0,934 mL/g olarak sabit alınır. Dengedeki iki faz için çözücünün kimyasal potansiyeli her iki uçta aynıdır; toplam bileşim bu aralıkta değişirken hesaplanan basınç bir plato verir. Bu sıkıştırılamaz model sabit referans hacmi yaklaşımını kullanır.

`ideal-osmotic.html` nötr ve moleküler olarak çözünmüş başka bir numune için seyreltik limitte RT/Mₙ eğimini ve seçilen derişimdeki ideal katkıyı hesaplar. Kullanıcı ölçülmüş Mₙ veya M𝑤 ve Đ girebilir. Yapı atölyesinden gelen sonlu molekülün mol kütlesi ancak kullanıcı bunu tek bir moleküler tür olarak seçerse hesaba girer.

### Deneysel karşılaştırma

`transfer.html` aynı χ = A + B/T + Cφ serbest enerji modelinin iki uyumunu karşılaştırır. Hacim kesrine bağlı χ için seyreltik ikinci virial katsayısı A₂ = (v̄²/Vₛ)(½ − A − B/T + C) olur. 73 bulutlanma sıcaklığı her iki uyumda ortaktır. Birleşik uyuma üç PS fraksiyonundan beş A₂ değeri eklenir. İki veri grubunun ortalama kare hataları kendi ölçeklerinde normalize edilir (1 K ve 10⁻⁵ mL mol/g²); ana karşılaştırmada göreli ağırlık λ=1, duyarlılık incelemesinde λ=0,25 ve 4 kullanılır. Ölçekler ve ağırlıklar modelleme tercihleridir, deneysel standart hatalar değildir.

| Dış sınama | Yalnız bulutlanma | Bulutlanma + ozmometri |
| --- | ---: | ---: |
| 20 ölçülmüş eşfaz çiftinde iki fazı veren koşul | 4/20 | 20/20 |
| Tüm 20 koşulda eşfaz genişliği MAE (hacim kesri) | 0,05157 | 0,02666 |
| İki modelin de eşfaz verdiği aynı 4 koşulda iki uç bileşimin RMSE’si | 0,04197 | 0,04547 |
| Bağımsız yayındaki 24 bulutlanma sıcaklığında RMSE | 0,774 K | 1,449 K |
| 8 membran ozmometrisi noktasında basınç MAE | 7,636 kPa | 10,063 kPa |

20 eşfaz çifti bir tek numuneye aittir (M𝑤 = 1,560 MDa, Mₙ = 1,515 MDa, toplam φ = 0,0324); 20 bağımsız polimer değildir. Tüm koşullarda deneyde iki faz bulunduğu için yanlış pozitif oranı ölçülmez. Tek faz öngörülen bir noktada genişlik sıfır alınır, var olmayan uç bileşimleri için sayısal değer atanmaz. Bulutlanma ve ozmometri, tek bir evrensel performans ölçütünde birleştirilmemiştir. Sonlu derişimdeki sekiz basınç ve üç dış A₂ ölçümü farklı numune ve yayınlardan gelir. Dış basınç ölçümleriyle uyum yapılmamıştır.

- [McIntyre ve ark. · kalibrasyon ozmometrisi](https://doi.org/10.1021/ja01523a003)
- [Nakata ve ark. · eşfaz eğrisi](https://doi.org/10.1103/PhysRevA.18.2683), [Wohlfarth derlemesi, 1978NAK kaydı](https://www.routledge.com/CRC-Handbook-of-Liquid-Liquid-Equilibrium-Data-of-Polymer-Solutions/Wohlfarth/p/book/9781032099675)
- [Striolo ve ark. · bağımsız A₂](https://doi.org/10.1021/ma000348e), [Nakamura ve ark. · bağımsız teta sıcaklığı](https://doi.org/10.1021/ma00017a029)
- [Krigbaum ve Geymer · bağımsız sonlu derişim basınçları](https://doi.org/10.1021/ja01517a022)

### Hesabı tekrarlama

Aşağıdaki komutlar Sites kaynağının ana dizininde çalışır. Düz GitHub Pages dağıtımı için dosya yollarındaki `dist/` ön ekini kaldırın. NumPy ve SciPy fizik modelini yeniden uyarlar; şekil için ayrıca Matplotlib gerekir.

```sh
python scripts/fit_transfer.py --cloud dist/ps-cloud-points.csv --osmotic dist/transfer-inputs/osmotic.json --coexist dist/transfer-inputs/coexist.json --out dist/transfer-results.js
python tests/test_transfer.py
node tests/test_osmotic.cjs
python scripts/plot_transfer.py dist/transfer-results.js dist/PS_Siklohekzan_Dis_Sinama.png
```

## Derişik çözeltide basınç ve çözücü aktivitesi

`dense.html` PS–siklohekzan için Mₙ = 72 kg/mol numunesinin basınç izotermlerini açar. Temel modelin A, B, C katsayıları önceki sürümdeki bulutlanma ve A₂ uyumundan alınır. İntegral χ'ye eklenen Dφ² terimi her veri bölmesinde yalnız eğitim basınçlarına uydurulur. 24 °C serisi, 34 °C serisi veya iki izotermin en yüksek derişimleri sırayla dışarıda bırakılır. Sekiz noktanın tamamına uydurulan katsayı ayrıca 44 °C'de sınanır. Bu seride iki kaynakta doğrulanan bir derişim ana sonuçtur. Kalan iki derişimin ikincil kaynağa dayandığı duyarlılık görünümünde belirtilir.

`solvent-activity.html` aynı fraksiyonda izotermal damıtma ile elde edilen, eğitimde kullanılmamış beş aktiviteyi inceler. Saf siklohekzan referanslı aₛ, −ln aₛ ve Δμₛ arasında geçiş yapılabilir. Bu karşılaştırmada ölçüm yöntemi ve bileşim aralığı birlikte değişir. Kaynak kayıtları ve belirsiz basamaklar `transfer-inputs/ACTIVITY.md` ile `transfer-inputs/PRESSURE44.md` dosyalarında açıklanır.

Aday modelin bağımsız θ bölgesi verisiyle çelişen negatif üçüncü virial katsayısı nedeniyle önceki uygulamanın varsayılanı korunur. A₃ ≥ 0 işaret kısıtı çoğu uydurmada A₃ = 0 sınırına gelir. Bu sınır da ölçülen pozitif büyüklüğün yeniden üretildiğini göstermez. İzoterm eğrileri denge serbest enerjisinden hesaplanır. Modelin iki fazlı bölgesinde ortak teğet platosu oluşur. Basınç noktalarının deneysel faz sayısı bilinmiyor.

```sh
node scripts/fit-dense.cjs
node tests/test-dense.cjs
node tests/test_activity.cjs
```

İsteğe bağlı Eφ³ duyarlılık hesabı `python scripts/probe-e.py` komutuyla yeniden üretilebilir (NumPy ve SciPy). Bu model, Dφ² sonuçları incelendikten sonra tasarlandı.

## Zincir boyutu ve örtüşme

- **chain.html** altı lineer PS fraksiyonunun siklohekzanda 34,5 °C ölçülmüş seyreltik Rg verisini kullanır. Bir çalışma içinde fraksiyon dışarıda bırakma sınaması, ayrı yayından düşük kütleli iki SAXS ölçümü ve indirilebilir ebeveyn-seri sınaması ayrı raporlanır. Model seçenekleri logaritmik uyum ve sabit θ üssüdür. İstenilen başka bir numunede ölçülmüş Mw ile aynı koşuldaki Rg verilirse geometrik örtüşme hesabı yapılır. Sıcaklık, numune ve mol kütlesi türü önceki uygulamalardan aktarılır.
- **chain-solvents.html** aynı altı fraksiyonun benzendeki 25 °C ve siklohekzandaki 34,5 °C ölçümlerini kimlik üzerinden eşler. Her iki ortamdaki Rg, Mw, geometrik c* ve seçilen derişimdeki c/c* gösterilir. Kaynaklar Miyaki, Einaga ve Fujita'nın Tablo I'i ile sayıları yeniden basan Fetters ve arkadaşlarının Tablo 7 ve 8'idir. Karşılaştırma çözücü ile birlikte sıcaklığın da değiştiği bu koşullara aittir.
- **chain-rotation.mp4** 64 adımlı tek bir θ Gauss konformasyonunun 8 saniyelik rijit dönüşüdür. Uzunluk ölçeği IK1500-1'in ölçülen kareli-ortalama Rg değeridir. Hesap, tek tek konformasyonların boyutunu sonradan bu değere eşitlemez. Tarayıcıdaki etkileşimli zincir görünümü farklı tohumlardan yeni bağımsız örnekler oluşturur. Zaman ölçeği veya kinetik parametre uydurulmaz.
- **solvent-activity.html** basınç uyumuna alınmayan iki ek ozmometri fraksiyonunun 13 damıtma kaydına aktarım sonuçlarını içerir. Bu numuneler saçılma fraksiyonlarıyla eşleştirilmez.

Kaynak izleri, ham saçılma kayıtları, sayısal kontroller ve birimler `chain-inputs/` altında sunulur. `node tests/test-chain.cjs` sayısal kontrolleri çalıştırır. `node scripts/chain-validation.cjs` doğrulama raporunu, `node scripts/check_mass_activity.cjs` önceki özellik uygulamasının damıtma sınamasını üretir. `python scripts/render_chain_rotation.py` MP4'ü yerel Node, Pillow ve FFmpeg ile yeniden oluşturur.

## Çözelti viskozitesi ve makine öğrenmesi

**viscosity.html** ölçülmüş intrinsik viskoziteyi [η] (mL/g), kapiler bağıl viskoziteyi ηᵣ ve koşula özgü saf çözücü korelasyonu mevcutsa çözeltinin dinamik viskozitesini (mPa·s) ayrı gösterir. PS’nin aynı altı fraksiyonu için 34,5 °C siklohekzan ve 25 °C benzen ölçümleri zincir boyutu uygulamasıyla bağlantılıdır. İki çözücü için ayrı Mark–Houwink uyumları, aynı çalışmanın ek üç fraksiyonunda ve farklı çalışmadan gelen altı tam sıcaklık eşleşmesinde sınanır. [η]c ≈ 1 operasyonel bir derişim ölçeğidir. Huggins katsayısı bulunmayan bir numunede seyreltik ikinci mertebe tahmini üretilmez. Genel bir yapıyla girişte kullanıcıya kendi ölçümünü aktarma olanağı sunulur.

**viscosity-transfer.html** aynı fraksiyonun ölçülen iki Rg değeri ve kaynak [η] ile hedef [η] tahminini hesaplar. Siklohekzandan benzene aktarımda yalnız Rg³ yaklaşımı altı fraksiyonda ortalama %34,96 mutlak sapma verir. Tüm bir ana numunenin fraksiyonları saklanıp kalan ana numuneden öğrenilen düzeltme kullanıldığında bu değer %3,86 olur. İki ana numune aynı kaynak serisindendir. Ölçülmüş hedef Rg, her iki tahminde de girdi olarak kullanılır.

**viscosity-ai.html** bir Huggins ortalama fonksiyonu ile u = c[η] üzerinde Gauss süreci düzeltmesini yan yana çalıştırır. Dört polistiren kapiler serisinde 20 ölçüm vardır. Her bölmede eğitim ve saklı gözlemler kimlikle ayrılır. PS–siklohekzan için bütün bir mol kütlesi dışarıda bırakıldığında, her iki numunenin [η] değeri bağımsız Einaga serisindeki kuvvet yasasından ekstrapole edilir. İki yönlü sınamada bağıl RMSE fiziksel modelde %0,604/%0,585, Gauss süreci eklendiğinde %0,411/%0,382 olur. Buna karşılık iki toluen serisinin en yüksek derişiminde AI ilavesi hatayı artırır. Dokuz çekirdek seçiminin sonradan yapılan duyarlılık karşılaştırması ayrı JSON dosyasındadır. Arayüzdeki tekil RMSE boyutsuz ηᵣ birimindedir. CSV'deki GP standart sapmaları sabit hiperparametreler altında gizli fonksiyona koşulludur.

**viscosity-formulation.html** tuzlu suda %5 mol hidrolizli HPAM'ın 24 °C'deki görünür viskozitesini polimer derişimi (kaynağın ppm birimi), NaCl (g/L) ve kayma hızı (s⁻¹) ile hesaplar. Akbari ve arkadaşlarının 20 deneyi 15 benzersiz durumu ve 9 reçeteyi kapsar. Her turda aynı derişim ve tuz miktarındaki bütün hızlar ve tekrarlar birlikte ayrılır. Test ölçütleri benzersiz durum başına eşit ağırlıkla fiziksel mPa·s cinsinden hesaplanır. Kullanıcıya gösterilen son modelin çalışma alanı ölçüm durumlarının üç boyutlu dışbükey zarfıyla sınırlandırılır. Kimyalar arası viskozite aktarımı için bağımsız polimer ailelerinden ölçüm gerekir.

Birincil kaynaklar ve birimler `viscosity-inputs/` altındadır. Fizik, AI ve aktarım testleri ile grafikler aşağıdaki komutlarla yeniden üretilebilir. Sitedeki `dist` dizini GitHub Pages ve çevrimdışı paketin ana dizinine karşılık gelir. Bu düzenlerde komutlardaki `dist` yerine `.` kullanılır.

```sh
node tests/test-viscosity.cjs
node tests/test-viscosity-gp.cjs
node tests/test-formulation.cjs
node tests/test-viscosity-transfer.cjs
node scripts/viscosity-validation.cjs
node scripts/viscosity-gp-validation.cjs
node scripts/viscosity-gp-sensitivity.cjs
python scripts/train-formulation.py dist/viscosity-inputs/hpam-akbari.csv dist/formulation-model.js
python scripts/plot-viscosity.py --dist dist --output dist/PS_Viskozite_Sinamasi.png
python scripts/plot-viscosity-ai.py --dist dist --output dist/PS_HPAM_AI_Sinamasi.png
```

Dokuz HPAM formülasyonu sırayla ayrıldığında, koşul başına RMSE kuvvet yasasında 1,781 mPa·s, ikinci derece Ridge modelinde 0,855 mPa·s ve aynı özniteliklerle cezasız ikinci derece modelde 0,685 mPa·s olur. Üç katsayı kümesi ve her katın hatası arayüzde karşılaştırılır.
