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
