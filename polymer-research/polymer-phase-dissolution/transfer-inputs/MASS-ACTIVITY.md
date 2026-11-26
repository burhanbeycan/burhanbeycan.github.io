# Farklı PS fraksiyonlarında siklohekzan aktivitesi

Krigbaum ve Geymer, *J. Am. Chem. Soc.* **81** (1959), 1859–1868. Basılı sayfa 1864, Tablo III. [Birincil yayın](https://doi.org/10.1021/ja01517a022) · [Tabloları içeren makale](https://electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20the%20American%20Chemical%20Society%20US/1959%20%20%28vol%20081%29/08%20%20%281777-2030%29/1859-1868.pdf)

Tablo III polistiren–siklohekzan numunesi ile kalibre edilmiş trifenilmetan–siklohekzan referansı arasındaki izotermal damıtma dengesini verir. Fraksiyon I ve III'ün tablodaki tüm satırları, kaynakta yazılı hacim kesri, sıcaklığı ve sayısal basamakları ile CSV'ye aktarıldı. Dördüncü sütundaki referans trifenilmetan bileşimi polimer hacim kesrinden ayrıdır. Hesapta beşinci sütundaki −ΔF̄₁/(RT) = −ln(aₛ) ve aₛ = exp[−(−ln aₛ)] kullanılır. Yazarların trifenilmetan referansı için verdikleri sıcaklığa bağlı kalibrasyon bağıntısından bağımsız olarak hesaplanan −ln(aₛ), aktarılan 13 noktada tablodan en çok 0,00057 ayrılır. Bu fark sayısal aktarım için bir tutarlılık denetimidir. Aktarım arşivdeki sayısal tabloya dayanır; satır düzeyinde bir deneysel standart hata yayımlanmamıştır.

Yazarların *Materials* bölümünde ve Tablo I'de bildirdikleri ozmometrik sayıca ortalama mol kütleleri Fraksiyon I için 25.100, II için 72.000 ve III için 440.000 g/mol'dür. Manometriyi anlatan Tablo IV'ün erişilebilen sayısal başlığında I için 25.900 g/mol okunur. Ana hesapta numune açıklamasını esas aldık ve Fraksiyon I için 25.900 g/mol ile hesaplanan sonuçları ayrıca yayımladık. Üç modelde ve dört noktada iki mol kütlesi arasında görülen en büyük öngörü farkı aktivite biriminde 0,000081'dir. Kullanılan sayı ortalama mol kütlesi Mₙ'dir. Bu ozmometri numuneleri için saçılmadan alınmış M𝓌 ve Rg bu veri setinde bulunmuyor.

Temel modelin A, B, C katsayıları önceki bulutlanma ve A₂ uyumundan sabittir. Ek Dφ² teriminin katsayısı, yalnız Fraksiyon II'nin 24 ve 34 °C'deki sekiz ozmotik basınç kaydından alınır. Aynı kaynaktaki 34 ve 44 °C damıtma verilerinin basınç uyumu dışında bırakılması, gözlem yöntemi ve numuneler arasında yapılan geriye dönük bir aktarım kontrolüdür. Her fraksiyonda yalnız Mₙ değiştirilir; bütün model katsayıları ve referans hacimleri sabittir. Denge hesabı faz ayrımı varsa ortak teğeti kullanır. Fraksiyon I ve III farklı bileşimlerde ölçüldüğünden fraksiyonlar arası hata karşılaştırması mol kütlesinin tek başına nedensel etkisini ölçmez. Bu fraksiyonlar zincir boyutu çalışmasındaki saçılma numunelerinden ayrıdır.

| Dış fraksiyon | Damıtma satırı | aₛ ortalama mutlak hata, önceki model | D serbest | A₃ ≥ 0 kısıtlı |
| --- | ---: | ---: | ---: | ---: |
| I · 25,1 kg/mol | 4 | 0,01160 | 0,00927 | 0,01058 |
| III · 440 kg/mol | 9 | 0,00823 | 0,00674 | 0,00749 |

Serbest D her iki grupta da ortalama mutlak aktivite hatasını düşürür. Aynı D değerinin θ bölgesinde öngördüğü negatif üçüncü virial katsayısı, bağımsız ışık saçılması sonuçlarıyla çatışır. Dolayısıyla model araştırma karşılaştırmasında kalır. Birbirine yakın noktalarda ölçülen aktivite 1'e çok yaklaşabildiğinden hesap dosyasında −ln aₛ hatası da yayımlanır. [Üçüncü virial verisi](https://doi.org/10.1021/ma00017a029)

`node scripts/check_mass_activity.cjs` komutu mevcut termodinamik motoru ve kayıtlı basınç uyumunu yükler, ham CSV'den model tahminleri ile hata ölçütlerini `transfer-inputs/mass-activity-results.json` dosyasında yeniden üretir.
