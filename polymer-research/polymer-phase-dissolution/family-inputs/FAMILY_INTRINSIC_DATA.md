# Teta koşulundaki aile aktarım sınamasının ölçümleri

Her satır, belirli bir polimerin bir mol kütlesinde, belirtilen çözücü ve sıcaklıkta ölçülmüş intrinsik viskozitesini içerir. `eta_mL_g`, sonsuz seyreltmeye ekstrapole edilen [η] değeridir. 1 dL/g = 100 mL/g. Mol kütlesi birimi g/mol'dür. Her çalışmanın kendi teta veya tetaya yakın ölçüm koşulu korunmuştur. Aynı aileyi dışarıda bırakan sınama, çoğu durumda o aileye özgü çözücü ve çalışma koşullarını da dışarıda bırakır. Ayrı bir sıcaklık ya da çözücü etkisi bu düzenekten belirlenemez.

| Polimer | Aile | Satır | Ortam | Ölçüm ve seçim |
|---|---|---:|---|---|
| Ağırlıkla 1,4-polibütadien (PBD) | Polidien | 9 | 1,4-dioksan, 26,5 °C | Fetters Tablo 1, özgün kaynak 23, hem M hem [η] bulunan satırlar. |
| Ağırlıkla 1,4-poliizopren (PI) | Polidien | 19 | 1,4-dioksan, yaklaşık 34,5 °C | Fetters Tablo 3, aynı özgün kaynak ve sütun ölçütü. Sıcaklık derlemenin seri için verdiği yaklaşık değerdir. |
| Poliizobütilen (PIB) | Poliolefin | 7 | İzoamil izovalerat, 25 °C | Abe Tablo I–II, 10.000 g/mol üstündeki, Đ ≤ 1,10 bildirilen numuneler. |
| Polistiren (PS) | Vinilaromatik | 12 | Siklohekzan, 34,5 °C | Einaga Tablo I, altı yüksek mol kütleli fraksiyon ve F serisinden üç standart. Baysal Tablo 1–2, Đ ≤ 1,10 bildirilen üç standart. |
| Poli(α-metilstiren) (PAMS) | Vinilaromatik | 9 | Siklohekzan, 37 °C | Cowie Tablo 1, dokuz numune. Özgün çalışmada ikinci virial katsayısından θ = 37,0 ± 0,5 °C bulunmuştur. |

`mass_basis=Mw` özgün yayındaki ağırlık ortalamalı mol kütlesidir. Polidien derlemesi mutlak mol kütlesi tayini ölçütü uygular ve tabloda M yazar. Bu nedenle bu satırlar `M_abs` olarak işaretlenir. Derlemenin sayısal kütle sütunundan oluşturulan `M0.108` gibi polidien numune anahtarları makaledeki fraksiyon adı değildir. `sample` ve `source` birlikte numuneyi tanımlar. Boş bir Đ hücresi, o numune için sayısal bir değer aktarılmadığını gösterir. PAMS çalışmasında seçili örneklerde Đ < 1,10 bildirildiğinden dokuz satıra ayrı ayrı 1,10 yazılmamıştır.

`repeat_smiles`, temel tekrar biriminin zincirdeki bağlantısını verir. Polidien numunelerinin küçük 1,2/3,4 katkıları, cis/trans dağılımları, uç grupları ve stereokimyası bu şemada çözülmez. Aile sınıfları ölçülen [η] değerinden bağımsız olarak polimer kimyasına göre tanımlanmıştır. Bu kohort beş hidrokarbon polimerini ve yalnızca üç geniş aileyi temsil eder. Aynı ana polimerin fraksiyonları ve aynı yayındaki ölçümler ortak kökene sahip olabilir.

`table` ve `ref`, Fetters derlemesindeki tablo ve özgün kaynak numarasıdır. Karşılaştırma amacıyla verilen derleme numarasının aynı çalışmanın yalnızca bir kısmını içerebildiğine dikkat edilmelidir. Örneğin PAMS numune 2'nin değeri ve küçük PIB numuneleri özgün tablodan alınmıştır. Kaynak numarasının boş olduğu Baysal satırları doğrudan özgün tablolara dayanır. PS ölçümleri önceki PS modülünden gelir ve burada her biri bir kez kullanılır. Ölçülen mol kütlesi ve [η] dışında yapı faktörü, benzetilmiş viskozite veya Mark–Houwink bağıntısından üretilen nokta eklenmemiştir.

## Kaynaklar

- İki polidien serisinin deneyi ve tablo derlemesi. [Hadjichristidis, Zhongde, Fetters ve Roovers (1982)](https://doi.org/10.1002/pol.1982.180200415) · [Fetters ve ark. (1994), Tablo 1 ve 3](https://www.nist.gov/system/files/documents/srd/jpcrd479.pdf)
- PIB'nin mol kütlesi, dispersitesi ve izoamil izovalerattaki intrinsik viskozitesi. [Abe, Einaga ve Yamakawa (1991), Tablo I–II](https://doi.org/10.1021/ma00015a026)
- PAMS'nin mol kütlesi ve 37 °C'deki çözücüye özgü viskozitesi. [Cowie, Bywater ve Worsfold (1967), Tablo 1](https://doi.org/10.1016/0032-3861(67)90014-6)
- PS'nin siklohekzandaki intrinsik viskozitesi. [Einaga, Miyaki ve Fujita (1979), Tablo I](https://doi.org/10.1002/pol.1979.180171205) · [Baysal ve Uyanık (1992), Tablo 1–2](https://doi.org/10.1016/0032-3861(92)90695-S)
