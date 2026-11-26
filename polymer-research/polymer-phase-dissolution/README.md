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
