# PEtOx uç modelleri için karşılaştırma

Her sistem H/H veya Me/OH uçları taşıyan bir DP 3 PEtOx zinciri ve 500 TIP3P sudan oluşur. Karşılaştırma kodu sıcaklık, basınç, NPT entegrasyon parametreleri, su kuvvet alanı dosya özeti, topoloji ve yörünge kimlikleri ile üretim zamanlarını doğrular. Her yeni üretim koşusuna çözücü ve zincirin kendi önceki son konum ve kutu koordinatları taşınır. Başlangıç hızları yeniden çekilir; yeniden ısıllaştırma üretimden ayrıdır. İki uç modelinde aynı tohum etiketleri kullanılmıştır. Farklı sayıda atom taşıyan sistemlerin kimyasal karşılaştırması üçer koşunun ortalaması üzerindendir. H/H modelinde önceki ve uzatılmış koşuların soyu kaynak yörüngesinin SHA-256 özetiyle ayrıca eşleştirilir.

## Koordinatlardan ölçülenler

1. Zincir, bağ grafındaki her adım için ortorombik kutuda en yakın görüntü kuralı uygulanarak bütünleştirilir. Tüm zincirin kütle ağırlıklı jirasyon yarıçapı hidrojenleri içerir. Ek uç atomlarının bu ölçüye kütle katkısı vardır. Ortak ağır atomlar üzerinde ikinci bir jirasyon yarıçapı hesaplanır. H/H yapısında C15N3O3 ağır atomlarının tamamı alınır. Me/OH yapısında α ucundaki metil karbonu ve ω ucundaki hidroksil oksijeni çıkarılınca aynı C15N3O3 çekirdeği kalır. Zincirin geri kalanında son konformasyon üzerindeki uç etkileri korunur.

2. Üç amid karbonili kimyasal bağlanmaya göre atanır. α ucundaki amid azotu H/H zincirinde H, Me/OH zincirinde metil taşır. Diğer iki amid azotu bu azottan sırasıyla üç ve altı bağ uzaktadır. Her oksijene 0,35 nm mesafe içinde bulunan su oksijenleri minimum görüntü kuralıyla sayılır. Bir su birden fazla karbonile komşu olabilir. Her karbonile suyun verdiği hidrojen bağı için ek ölçüt su O–H···karbonil O açısının en az 150° olmasıdır. Sayım hidrojen verici–H–alıcı üçlüleri üzerindendir. Terminal N–H veya O–H grubunun suya verdiği bağlar ayrı bir sütundadır. Su → terminal hidroksil oksijeni yalnız bu grup mevcut olduğunda hesaplanır.

3. Her koşu için gerçek zaman damgaları, üçer bağımsız koşudan elde edilen ortalamalar, dört ardışık bloğun ortalamaları ve ikinci yarı eksi ilk yarı özeti kaydedilir. Seriler eşit fiziksel süreye ve örnekleme aralığına sahiptir. Örneklem standart sapması üç koşu ortalamasından hesaplanır. Eski H/H koşularının karşılaştırması için kendi özgün üretim aralığı ve kayıt sıklığı kullanılır. H/H modelinin erken ve geç yeni pencereleri ayrı tutulur; çalışma süreleri birleştirilmez.

4. Her gözlenebilir ve koşu için örnek otokorelasyon dizisinden başlangıçtaki pozitif, monotonlaştırılmış ardışık çift toplamları ile istatistiksel verimsizlik **g** tahmin edilir. **g** değeri 1 ile kayıt sayısı arasında sınırlanır. Etkin örnek sayısı N/g, bütünleşik otokorelasyon zamanı g Δt/2 olarak tanımlanır. Sabit serilerde sonuç tanımsızdır. Yavaş ve ziyaret edilmemiş konformasyonlar bu ölçülerde ortaya çıkmayabilir.

5. Üç koşu yarıya bölünür. Birleştirilmiş sıralardan normal skorlar hesaplanarak, ayrıca ortak medyandan mutlak sapmalar katlanıp sıralanarak iki split R-hat tahmini elde edilir. Yorumda üstteki değer dikkate alınır. Tüm seri sabit olduğunda R-hat tanımsızdır. Bu uygulama MCMC literatüründeki çoklu koşu tanısını MD zaman serilerinde bir örnekleme göstergesi olarak kullanır. Küçük bir değer tek başına termodinamik dengeye erişildiğini kanıtlamaz. SS ve tek koşu ESS değerleri moleküler modelin sistematik hata payını kapsamaz.

## Koşullar ve fiziksel yorum

Aynı su sayısı her zincir için sabit bir zincir/su sayısı oranı tanımlar. Me/OH uçları kuru zincirin kütlesini artırdığı için kütlece polimer yüzdesi iki sistemde tam aynılık göstermez. NPT kutularının denge hacmi de dalgalanır. Me/OH modelinin kısmi yükleri bütün molekül için yeniden atanmıştır. α ve ω uçları bir arada değiştirilmiştir. Her ucun ayrı nedensel katkısı ve uçlar arası etkileşim için Me/H ile H/OH kontrolleri gerekir. Sonlu bir DP 3 oligomerinden sonsuz polimer zincirinin denge çözünürlüğüne doğrudan geçilemez.

## Kaynaklar

- Grossfield ve ark. MD örnekleme kalitesi, otokorelasyon ve belirsizlik hesabı. [Living Journal of Computational Molecular Science](https://doi.org/10.33011/livecoms.1.1.5067)
- Vehtari ve ark. Sıralara dayalı ve katlanmış split R-hat yöntemi. [Bayesian Analysis](https://doi.org/10.1214/20-BA1221)

## Su molekülü değişimi

Aynı karbonilin çevresindeki komşu sayısı yaklaşık sabitken yakındaki suların kimliği değişebilir. Su oksijeni kimlikleri yörünge boyunca korunur. 1, 2, 5 ve 10 ps gecikmeler için ilk karede belirli karbonile komşu olan her su ele alınır. Oranların paydası, karşılaştırılan gecikme için son karesi mevcut olan başlangıçlarda sayılan tüm (başlangıç karesi, su kimliği) komşuluklarıdır. İlk oran aynı suyun iki uç karede komşu olmasını gerektirir. İkinci oran aradaki kaydedilmiş her karede de komşuluk arar. Süreklilik yalnız kaydedilmiş kareler için gözlenmiştir. Koşu başına oranlar eşit ağırlıkla birleştirilir. Farklı zaman başlangıçları korelasyonludur. Bu analiz suyun izlenen yerde kalması hakkındadır; tekil bir hidrojen bağının kinetik ömrünü vermez.

Zincir boyutlarının dağılımını karşılaştırmak için aynı sayıda kare içeren iki koşunun gözlenen değerleri sıralanır. Eşleşen sıralı değerler arasındaki ortalama mutlak fark bir boyutlu ampirik Wasserstein-1 uzaklığıdır. Üç koşu çiftinden en büyüğü kayan dağılımları özgün birimde (örneğin Rg için nm) gösterir. Zaman noktalarına bağımsız gözlem varsayımıyla p değeri atanmaz.
