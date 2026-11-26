# θ çözücülerinde aile aktarımı

## Sınanan hedef

Girdi polimerin nominal tekrar birimi ile kaynakta bildirilen mol kütlesidir (g/mol). Hedef, o polimerin kaynakta θ koşulu olarak listelenen çözeltisinde intrinsik viskozitedir (mL/g). Her satır seyreltik limitteki [η] ölçümüdür. Bu sınama belirli bir çözücüde kendiliğinden θ sıcaklığı bulmaz. Derişime ve kayma hızına bağlı akış eğrisi başka bir hedeftir.

Üç dış kattan her birinde bir kimyasal iskelet grubunun bütün satırları saklanır. Polidien grubu 1,4-polibütadien ve 1,4-poliizopreni, vinil aromatik grup polistiren ve poli(α-metilstiren)i, üçüncü grup alifatik poliizobütileni içerir. Böylece yakın analoglar bir katın iki tarafına dağılmaz. Çözücü ve sıcaklık kimyaya göre değişir; dış sınama kimya ile deney koşulunun birlikte aktarımını ölçer. Kaynak kimlikleri aynı derlemede ilk yayına kadar izlenir. Derleme tek başına bağımsız dış laboratuvar sınaması oluşturmaz. Bildirilen M türü satırla birlikte korunur.

## Denklemler

İdeal uzun zincir için θ koşulunda Rg² ∝ M ve hidrodinamik hacim yaklaşımında [η] ∝ Rg³/M bağıntıları kütle üssünün 1/2 seçimine fiziksel dayanak sağlar. Sonlu zincirler, mol kütlesi dağılımı ve hidrodinamik etkileşimler sapmaya yol açabilir.

Temel model

    y = ln([η]/(1 mL/g)) − 0,5 ln(M/(100 000 g/mol))
    y_tahmin = b

Tanımlayıcılı model

    x = ln(j / (M₀ / (1 g/mol)))
    z = (x − eğitim_ortalaması) / eğitim_standart_sapması
    y_tahmin = b + βz

j nominal tekrarlayan birimin zincire eklediği iskelet bağı adımıdır. Omurgadaki çift bağ bir adım sayılır. 1,4-dienlerde j=4, diğer üç kimyada j=2 alınır. M₀ tekrar birimi formülünden C için 12,011 ve H için 1,008 g/mol kullanılarak hesaplanır. Nominal 1,4 bağlantısı küçük yan izomer paylarını ve cis/trans dağılımını ortalamaz. Freely jointed zincirde bağ sayısı yaklaşık jM/M₀ ile artar. Gerçek karakteristik oran, etkin bağ uzunluğu ve çözücüye bağlı yerel konformasyon da Rg’yi belirler. Öğrenilen tek katsayı bu fiziksel büyüklükleri ayrıştırmaz.

Her eğitim ailesinin ağırlığı eşittir; bir ailedeki kimyalar, bir kimyadaki birincil kaynaklar ve kaynak içindeki gözlemler kendi düzeylerinde eşit ağırlığa sahiptir. Toplam eğitim ağırlığı birdir. Amaç Σwᵢ(yᵢ−b−βzᵢ)² + λβ² olup sabit λ=1 ve serbest kesişim kullanılır. Bu değer, eğitim varyansı bire ölçeklendiğinde cezasız eğimi yarıya küçültür. Sınama sonuçları kullanılarak ayar seçilmez. Her iki model aynı eğitim satırları, fiziksel kütle üssü, yanıt ve ağırlıkları alır. Tanımlayıcı doğrudan nominal yapıdan gelir.

## Karşılaştırmayı okuma

Her test noktasına tekil deney kimliği, birincil kaynak kaydı, mol kütlesi türü, çözücü ve sıcaklık açıklaması eşlik eder. Hem kütle eksenindeki hem tanımlayıcı eksenindeki eğitim aralığı ayrıca kontrol edilir. Bu eksenlerde değer aralıklarının çakışması, saklı kimyanın eğitimde temsil edildiğini göstermez. Üç ailelik protokol için yeni ailelerde kapsaması ölçülmüş bir tahmin aralığı sunulmaz.

- MAPE her noktanın |tahmin−ölçüm|/ölçüm oranının yüzde cinsinden ortalamasıdır.
- Aile dengeli MAPE üç aile MAPE değerinin aritmetik ortalamasıdır.
- Aile dengeli bağıl RMSE'de önce her ailede bağıl hataların karesi ortalanır, sonra üç aile eşit ağırlıkla ortalanır ve karekök alınır.
- Birleşik puanlar her gözleme eşit ağırlık verir. Etkin kimyasal çeşitlilik beş kimya ve üç iskelet grubudur.

Birincil verilerin birimleri ve kökeni [Fetters ve ark. (1994)](https://doi.org/10.1063/1.555949) ile; modelin her katını ve sonuç satırlarını uygulamayla dağıtılan CSV/JSON çıktılarıyla karşılaştırabilirsiniz.

## Kayıt niteliği ve duyarlılık

28 polidien satırı özgün kaynakta mutlak mol kütlesi olarak bildirilir (M_abs). Diğer 28 değer M𝑤 olarak bildirilmiştir. M_abs ayrıca bir dağılım ortalamasına dönüştürülmedi. Yalnız açıkça M𝑤 olarak bildirilen satırlar tutulursa üç iskeletli değerlendirme gerçekleştirilemez. Polidienlerin 19 PI satırı için tek numuneye özgü sıcaklık derlemede yaklaşık olduğundan ölçüm arayüzünde bu niteliğiyle işaretlenir.

`theta-sensitivity.cjs` ana sonuçlardan sonra tanımlanmış iki duyarlılık incelemesi uygular. Birincisi yaklaşık θ sıcaklığına sahip PI satırlarını çıkarır. İkincisi üç ailenin ölçülen mol kütlesi aralıklarının kesişiminde kalır. Her iki yaklaşım da üç dış katı, tüm eğitim ağırlıklarını ve katsayıları altkümede yeniden oluşturur. Sonuçlar, her grubun örnek sayısı ile birlikte `family-inputs/theta-sensitivity.json` içinde yer alır.
