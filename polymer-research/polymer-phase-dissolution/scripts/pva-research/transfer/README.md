# Sulu çözeltilerde polimerler arası viskozite aktarımı

Dört malzeme için yayımlanmış 286 farklı deney koşulu karşılaştırıldı. Hedef belirli bir sıcaklık ve kayma hızındaki görünür dinamik viskozitedir (cP = mPa·s). Her malzeme için bir yayın bulunuyor. Bir kimyayı dışarıda bırakmak aynı zamanda yeni bir laboratuvara, numuneye ve çoğu zaman yeni bir koşul aralığına geçmeyi gerektiriyor.

## Ölçümler

| Kimya | Benzersiz koşul / izotermal reçete | Özgün deney ve birim | Kaynak |
|---|---:|---|---|
| PVA | 244 / 16 | Merck Peru 141351, suda ağırlıkça %4–10, 20–35 °C. Deneysel τ (Pa) ve γ̇ (s⁻¹) çiftlerinden 1000τ/γ̇. Sayısal DH ve doğrulanmış Mₙ eksik. | Carrasco-Venegas vd., [Yöntem ve Tablo 2](https://doi.org/10.3390/polym15071743). CC BY 4.0. |
| HPAM | 15 / 9 | %5 mol hidroliz, sulu NaCl, 24 °C. Derişim makaledeki ppm biriminde. Yayımlanan log₁₀η etiketinden fiziksel η'ya dönüş. Aynı orta koşulun altı tekrarı cP cinsinden ortalama alınarak birleştirildi. | Akbari vd., [Tablo 1 ve 2](https://doi.org/10.1007/s13202-017-0375-6). CC BY 4.0. Girdi dosyasında 20 deney. |
| Karboksimetilselüloz (CMC) | 18 / 3 | Suda 4 g/100 mL, DS 0,43/0,88/1,23, 25 ± 1 °C. 170–1.000 s⁻¹. Tabloda doğrudan cP. Holoselüloz karışımları hariç. | Laribi vd., [Yöntem ve Tablo 4](https://doi.org/10.35812/CelluloseChemTechnol.2020.54.24), [PDF](https://www.cellulosechemtechnol.ro/pdf/CCT3-4%282020%29/p.225-236.pdf). |
| Sodyum aljinat (ALG) | 9 / 9 | Kuru katı kütlesi ve su kütlesinden hazırlanmış %1–3 çözelti. 16 h hidratasyon (20 °C). 25, 40 ve 55 °C'de 100 s⁻¹ için bildirilen görünür viskozite Pa·s × 1.000. Cross modeli çıktıları dahil edilmedi. | [Birincil makale](https://doi.org/10.1111/j.1365-2621.2009.02088.x) · [Tablo 1](https://academic.oup.com/view-large/495475116). |

Her satır özgün gözlemin DOI'sini ve tablosunu taşır. Sıcaklık kaydı kaynağın nominal değeridir. Kütle, kütle/hacim ve temeli belirtilmeyen ppm derişimleri ayrı özniteliklerde tutuldu. NaCl girdisi bildirilen ilave tuzdur. Toplam iyonik kuvvet ayrıca ölçülmemiştir. near_298.csv içindeki filtre **297,15–298,15 K** nominal sıcaklıkları içerir. Hedef tam 298,00 K'ye dönüştürülmemiştir.

## Model ve doğrulama

1. Dış bölme ya tek bir malzemenin bütün gözlemlerini ayırır (material, dört kat) ya da glikozidik polimerlerin ikisini birlikte ayırır (group, üç kat: PVA, HPAM, CMC+ALG). 
2. Koşul modeli log viskoziteyi log kayma hızı, 1/T, özgün derişim kanalları, derişim temeli ve ilave tuzla ilişkilendiren Ridge modelidir. Kimya modeli glikozidik omurga, amid içeren ana motif ve nominal tekrar birimi başına potansiyel karboksil grubundan bir benzerlik çekirdeği oluşturur. İyonlaşma ve gerçek yük ölçümü için ayrıca veri gerekir.
3. Önceden eğitilmiş polyBERT'in 600 bileşenli dondurulmuş çıktısıyla ek bir benzerlik çekirdeği hesaplanır. PVA ve HPAM için ana motif, CMC için karboksimetil ikamesi ayrıca girilen selüloz omurgası kullanılır. Aljinatın M:G oranı belirsiz olduğundan ana karşılaştırmada iki epimerin eşit ağırlıklı ortalaması tanımlayıcıdır. M ve G uç senaryoları her satırda ayrıca yer alır. [polyBERT](https://doi.org/10.1038/s41467-023-39868-6) · [Kullanılan ağırlık sürümü](https://huggingface.co/xushijie/polyBERT/tree/e7dce434fb3eff37905dc114008660e5479ca9a8).
4. Kimya çekirdeği bir ofset ve koşula bağlı eğim düzeltmesi sağlar. Gaussian çekirdek uzunluk ölçeği her dış katın yalnızca eğitimdeki benzersiz tanımlayıcı uzaklıklarıyla belirlenir. Sabit Ridge cezası 0,1. Eğitimde önce her malzeme, onun içinde her izotermal reçete, onun içinde her koşul eşit ağırlıktadır. Üç sabit ceza değeri üzerindeki duyarlılık ayrıca yayımlanır.
5. Tek-referans modu, her yeni **izotermal reçetenin** log kayma hızı aralığının ortasına en yakın tek ölçümü seçer ve tahmini o değere çarpanla bağlar. Yalnız geri kalan gerçek ölçümler puanlanır. Aynı test noktalarındaki referanssız sonuç zero_reference_on_anchor_test grubundadır. Sadece bir kayma hızı bulunan aljinat reçetelerinde bu ek değerlendirme uygulanamaz.

## Elde edilen sonuçlar

Üç ana kimya grubuna eşit ağırlık verilmiştir. Her grubun hatası kendi noktalarında hesaplandı. N = 286 nokta. log-MAE, |log₁₀(η̂/η)| değerinin ortalamasıdır.

| Model | Yeni grupta ortalama log-MAE | Yeni grupta ±%10 içindeki nokta oranı |
|---|---:|---:|
| Koşul modeli | 1,351 | %0 |
| Açık kimyasal tanımlayıcılar + koşullar | 1,068 | %0 |
| Dondurulmuş polyBERT + koşullar | 1,167 | %0 |
| Kimyasal tanımlayıcılar + polyBERT + koşullar | 1,001 | %0,27 |

Birleşik modelde bütün dış noktalarda ±%10 içinde kalan gözlem 2/286 (%0,70). Grup eşit ağırlıklı oran %0,27. PVA'yı tamamen ayırınca medyan mutlak yüzde hata %237,1 çıktı. CMC ve ALG birlikte ayrılınca medyan %97,3 oldu. Önceden eğitilmiş temsil, kimya modeline göre ortalama log hatasını kısmen düşürdü. Ceza 0,03/0,1/0,3 seçeneklerinde birleşik modelin grup düzeyinde ortalama ±%10 oranı %0–0,96 arasındaydı. Bu dış test rastgele bir numune örneklemesi değildir. Yeni polimer başarısı için geniş ve örtüşen koşullar, bağımsız laboratuvar ve tanımlı kütle serileri gerekli.

Geniş grup bölmesinde eşlenmiş 249 referans dışı nokta için, birleşik modelin ortalama log hatası referans ölçümsüz 1,051, her yeni izotermal eğri için bir gerçek referansla 0,0783 oldu. ±%10 içindeki grup ortalaması sırasıyla %0,29 ve %31,58. Bu, ölçümle koşullandırılmış eğri şekli denetimidir.

## Tekrar üretme

Aşağıdaki komutlar proje kökünden çalışır. Hesaplanmış 286 satırlı CSV ve donmuş gömmeler dahil. NumPy, pandas ve SciPy yeterlidir.

    python transfer/prepare_observations.py
    python transfer/benchmark.py
    python transfer/sensitivity.py
    python transfer/test_transfer.py

heldout_predictions.csv tüm dış tahminleri, varsa referans noktasının kimliğini ve aljinat temsil uçlarını içerir. by_fold.csv, by_curve.csv, macro_scores.csv ve near_298.csv dış grup, deney eğrisi, ortalama ve dar sıcaklık aralığı sonuçlarını verir. regularisation_sensitivity.csv aynı sabit katlarla alternatif Ridge katsayılarını verir.

## Kullanıcının yapısını araştırma modunda değerlendirme

Kayıtlı motif veya gerçek kodlayıcıyla oluşturulan yeni bir tekrar birimi kullanılabilir. Transfer modeli sulu çözelti için yapı ve reçeteyi kullanır. Eğitimde numunelerin Mₙ değerleri ortak biçimde ölçülmediği için bu modelde bir Mₙ katsayısı kestirilmemiştir. --Mn yalnızca istek üst verisidir. Çıktıda conditions_on_Mn: false belirir. Kütleye özgü PVA operatörü ayrı doğrulanmalıdır.

    python transfer/infer.py --material CMC --potential-carboxy-sites 0.88 --concentration 4 --basis g_per_100mL --T 298.15 --rate 300
    python descriptors/encode_polybert.py --model-dir /tmp/polybert --download --repeat '[*]CCO[*]' --out /tmp/peo-query.npz
    python transfer/infer.py --material userPEO --query-embedding /tmp/peo-query.npz --backbone other --amide-bearing-parent 0 --potential-carboxy-sites 0 --concentration 4 --basis wt_percent --T 298 --rate 10

Tek referans varsa --reference-rate ve --reference-cP birlikte girilir. Öngörü destek bilgisi tek tek deneysel değişkenlerin sınırlarını bildirir. Çok değişkenli ortak kapsama veya dış dağılımda deneysel olarak kalibre edilmiş tahmin aralığı sunmaz.
