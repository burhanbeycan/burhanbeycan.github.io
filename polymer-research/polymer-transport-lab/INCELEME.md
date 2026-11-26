# PVA filminde su difüzyonu

Neme bağlı bir difüzyon katsayısı su alımını hangi koşullarda sabit katsayıdan daha iyi açıklar?

## Asetat dizisi ve suyun zincir çevresindeki dağılımı

**Ölçüm.** Nakajima ve Gocho, düşük bağıl nemde homojen ve heterojen asetilasyonla hazırlanan örneklerin benzer miktarda su aldığını, yüksek nemde ise homojen hazırlanan serinin daha fazla su aldığını bildiriyor. Hazırlama biçimi ile kristal yapı birlikte değiştiğinden bunu tek başına dizi etkisi olarak atayamayız. İncelenen üç birincil kaynakta aynı numune için tam NMR, referans MD ve denge su alımı birleşimi doğrulanamadı. [Deneysel çalışma](https://doi.org/10.1246/nikkashi.1978.1431)

**Uygulamada dene.** “PVA · withhold complete grades (synthetic)” örneğinde her PVA sınıfı tamamen dışarıda bırakılır. Aynı 12 yapay sınıfta bileşim/koşul, ek NMR diyadı ve ek MD hidratasyonu sırasıyla 0,05490, 0,03863 ve 0,00465 g/g ortalama mutlak hata verir. Bu veriler bilinen bir formülden üretilmiştir. Gerçek tahmin için aynı numuneye bağlı iki CSV girilebilir. Kendi yörüngen için “PVA · oxygen coordination across a periodic boundary” seçeneği hidroksil ve asetat karbonili başına su komşularını sayar. [CSV şemaları ve değerlendirme](docs/PVA_SEQUENCE.md)

## Uygulamada incele

Denge izotermini ve kinetik ölçümleri ayrı gir. Bir veya iki açık yüz seç. Deneysel HDI/PVA örneğinde iki uç bağıl nemden ara nemi hesapla. Kalınlık sınamasında bir kalınlığa ait bütün eğriler eğitim dışına çıkar.

## Örnek sonuç

25 °C’de ve bağıl nem %61 iken yayımlanan görünür D, 1,40 × 10⁻¹¹ m²/s. Diğer iki nemden türetilen derişime bağlı model 1,51 × 10⁻¹¹, sabit model 1,67 × 10⁻¹¹ m²/s verdi.

![%61 bağıl nemde PVA için tek ölçüme karşı görünür difüzyon katsayısı tahminleri](docs/assets/benchmark-comparison.svg)

## Kimyasal yorum

Bu karşılaştırma bir ara nem noktasına aittir. Kalınlıklar arası örnek sentetiktir. Çözelti dökümünde kalınlığın zamanla değiştiği ve dış kütle aktarımının sınırlayıcı olabildiği durumlar ayrı bir kuruma problemi gerektirir.

[Bilimsel yöntem ve veri](docs/PVA_WATER.md) · [Proje ve kullanım](README.md) · [Tek proje kaynak arşivi](../downloads/polymer-transport-lab-0.4.0.zip) · [Tüm projeler](../README.md)
