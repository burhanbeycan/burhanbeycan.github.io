# Polistiren çözeltilerinde kapiler viskozimetri serileri

CSV’de her satır bir numune için tek derişimdeki kapiler deney sonucudur. İki toluen kaynağında akış süreleri, siklohekzan kaynağında hesaplanan dinamik çözelti viskozitesi bildirilir. `eta_rel` kaynağın yayımladığı boyutsuz orandır. Boş hücreler bildirilmeyen büyüklükleri gösterir.

| Seri | Numune ve koşullar | Ölçümlerin kaynağı |
|---|---|---|
| Copenhaver vd. | PS granül, tedarikçiye göre Mₙ = 170.000 ve Mw = 350.000 g/mol. Toluen, 30 ± 1 °C. | Beş ekip, her ekip için ayrı bir Ubbelöhde viskozimetresi ve üç tekrarlı ölçüm. Kaynaktaki gerçek derişim g/dL'den g/L'ye 10 ile çarpılarak çevrildi. İkinci ekibin sonucu yayın yazarları tarafından tartışılıyor. Tüm satırlar korundu. |
| PHYWE P3010601 | PS, toluen, 25,0 °C. Mₙ ve Mw bildirilmemiştir. | Eğitim amaçlı cihaz kılavuzunda her süre dört ölçümün ortalamasıdır. 1 mg/mL = 1 g/L. Kılavuzun aynı seriden ekstrapole ettiği [η] = 0,0644 L/g ve bundan türettiği Mᵥ ≈ 133.500 g/mol bağımsız numune değerleri olarak kullanılmaz. |

İki toluen kaynağı bağıl viskoziteyi kapiler akış süresi oranından tahmin eder. Aynı geometrik sabit ve ihmal edilebilir kinetik enerji düzeltmesiyle dinamik viskozitelerin oranı ρt/(ρ₀t₀) olur. Burada kaynakların seyreltik çözeltide kullandığı ρ ≈ ρ₀ varsayımı geçerlidir. Yayımlanmış üç ondalıklı ηᵣ değerleri kullanılır. PHYWE'nin indirgenmiş ve inherent viskozite sütunları yerine bu iki büyüklük doğrudan ηᵣ ve ölçülen derişimden hesaplanır.

İki toluen serisi farklı polistiren numunelerine, sıcaklıklara ve ölçüm protokollerine aittir. Her serinin derişimlere göre ayrılmış denetimi kendi numunesine aittir. Yüksek derişimde c[η] büyüdükçe Huggins veya Kraemer düşük derişim açılımlarının ekstrapolasyon hatası ölçülebilir. Bu tablolarda farklı çözücülerde aktarım ya da kayma hızı taraması bulunmaz.

## Siklohekzanda 34,5 °C

Luh’un Tablo 8’inin ikinci sayfasında yer alan iki tam derişim dizisi aktarıldı. Tablo 1 polistiren numunelerinin Mw ve Mw/Mn değerlerini veriyor. Kalibre edilmiş Cannon–Ubbelohde viskozimetrelerinde her derişim için akış süresi beş kez ölçülmüş. Dinamik viskozite hesabında kapiler kalibrasyon, kinetik enerji düzeltmesi ve çözelti yoğunluğu kullanılmış. Tablo 8 tek tek süreleri ve bunların standart sapmasını vermiyor. Ölçümlerde banyo sıcaklığının toleransı ±0,02 °C. Kütle derişimi, çözelti hacmi başına polimer kütlesidir. [Luh (1967), s. 14–22, 53–57](https://scholarsmine.mst.edu/masters_theses/6870/)

| Mw (g/mol) | Đ | Ölçülen c (g/L) | Ölçülen ηçözelti (mPa·s) | Tüm noktalardan hesaplanan [η] (mL/g) |
|---:|---:|---:|---:|---:|
| 633.000 | 1,06 | 3,371–9,956 | 0,9613–1,4702 | 66,22 |
| 1.770.000 | 1,08 | 1,157–4,704 | 0,8730–1,2907 | 110,89 |

Veri aktarımında Tablo 8’deki c değerleri g/dL’den g/L’ye 10 ile çarpıldı. Dinamik viskozitede 1 cP = 1 mPa·s. Her satırdaki çözelti viskozitesini aynı satırdaki bağıl viskoziteye bölmek yaklaşık 0,7673 mPa·s verir. Bu işlem aynı kapiler ölçümünden türeyen iki sütun arasındaki aritmetik tutarlılığı sınar. Beşer noktayla Huggins ekstrapolasyonları Tablo 5A’nın sırasıyla 66,2 ve 110,9 mL/g değerlerini basılı basamaklar içinde geri verir. Bu iki sayı aynı çalışmanın kontrolleridir. Çözeltinin bir kayma hızı taraması veya bağımsız tekrarlanabilirlik belirsizliği bu tablodan türetilemez. [Luh (1967), Tablo 5A, Tablo 8 ve Ek 2](https://scholarsmine.mst.edu/masters_theses/6870/)

## Kaynaklar

- Katie Copenhaver, Samuel Pennell, Samyak Jain, Paul S. Russo ve Prateek Verma (2024). *Classic Ubbelöhde Intrinsic Viscosity Laboratory Exercise Made Simple and Fast*. Journal of Chemical Education 101, 4952–4958. Tablo 1, ölçüm protokolü ve Şekil 3 tartışması. Makale CC BY 4.0 lisanslıdır. Bu kopyada tablo sütunları seçildi, derişim birimi dönüştürüldü ve satırlar sınıflandırıldı. [DOI](https://doi.org/10.1021/acs.jchemed.4c00439) · [Lisans](https://creativecommons.org/licenses/by/4.0/)
- PHYWE (P3010601). *Determination of the molecular weight of a polymer from intrinsic viscosity measurements*. Tablo 1, Şekil 2 ve kaynakta verilen ölçüm sıcaklığı. [Deney kılavuzu](https://spegroup.ru/upload/wikifiles%202/Experiment_guide_P3010601.pdf)

- Hwa Luh (1967). *The effect of concentration on the viscosity of dilute and moderately concentrated polymer solutions*, University of Missouri at Rolla, yüksek lisans tezi. Tablolar 1, 5A ve 8. [Üniversite arşivi](https://scholarsmine.mst.edu/masters_theses/6870/)
