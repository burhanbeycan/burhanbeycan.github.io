# Polistirenin intrinsik viskozitesi ve ölçülmüş zincir boyutu

Einaga, Miyaki ve Fujita, aynı altı lineer polistiren fraksiyonunun sıfır kayma hızındaki limit viskozite sayısını benzende 25 °C ve siklohekzanda 34,5 °C’de ölçtü. Çalışmada önceki saçılma deneylerinde kullanılan IK1500 ve BK2500 fraksiyonlarının adları açıkça verilir. Bu altı fraksiyon her çözücü için 8,76–56,8 milyon g/mol aralığında ayrı bir kalibrasyon sağlar. Aynı ana üründen ayrılan fraksiyonlar ortak köken taşır. [Einaga, Miyaki ve Fujita (1979), s. 2104 ve Tablo I, s. 2107](https://doi.org/10.1002/pol.1979.180171205)

## Hangi kayıtlar kullanılıyor?

| Seri | Her çözücüdeki numune sayısı | Kullanım |
|---|---:|---|
| IK1500-1, IK1500-2 ve BK2500-1 ila BK2500-4 | 6 | Kalibrasyon, her numunede aynı koşuldaki ışık saçılmasıyla eşleşme |
| F450, F288 ve F126 | 3 | Aynı yayından, her iki çözücüde ölçülmüş üç düşük mol kütleli F standardı. Kalibrasyonun alt sınırı altında ekstrapolasyon denetimi |
| Başka yayınlardan numuneler | Ayrı veri dosyasında | Bağımsız çalışma karşılaştırmaları için `viscosity-external-data.js` |

## Ölçüm ve birimler

- Viskozite yayınının Tablo I’inde iki çözücü için tek bir ağırlık ortalamalı mol kütlesi sütunu vardır. `Mw_g_mol` bu değerdir. Kalibrasyonda kullanılan bağımsız değişken Mw’dir. Viskozite ortalamalı mol kütlesi Mᵥ için ayrı bir ölçüm yapılması gerekir.
- Tablodaki intrinsik viskozite değerleri 100 cm³/g birimiyle basılmıştır. 1 dL/g = 100 mL/g. Hem ölçüm hem basılı ± değeri 100 ile çarpılmıştır. Örneğin IK1500-1 için benzende 5306 ± 50 mL/g, siklohekzanda 667 ± 7 mL/g bulunur. ± işaretine kaynakta istatistiksel güven düzeyi atanmamıştır.
- Huggins katsayısı `kH` boyutsuzdur. Yazarlar siklohekzanda en büyük üç mol kütlesi için katsayı yayımlamaz. İndirgenmiş viskozite–derişim eğrilerindeki yukarı yönlü kıvrım mol kütlesiyle belirginleşir. Bu üç kayıt eksik değer olarak saklanır. Bu tablodan her fraksiyon için sayısal bir derişim üst sınırı çıkarılmamıştır. [Birincil yayın, Şekil 3–4 ve Tablo I](https://doi.org/10.1002/pol.1979.180171205)

## Numunelerin saçılma ölçümüyle eşleştirilmesi

`Rg_nm`, aynı çözücü ve sıcaklıkta statik ışık saçılmasından elde edilen √⟨Rg²⟩z’dir. Anahtar, yayında tanımlanan fraksiyon adıdır. Saçılma yayını iki çözücüde ayrı Mw tayinleri bildirdiğinden `Rg_Mw_g_mol` ayrıca tutulur. Örneğin IK1500-1 için siklohekzandaki saçılma tayini 55,7 milyon g/mol, viskozite makalesinde kullanılan ortak değer 56,8 milyon g/mol’dür. F serisinden numunelere bir Rg atanmamıştır. [Miyaki, Einaga ve Fujita (1978), Tablo I, s. 1183](https://doi.org/10.1021/ma60066a022) · [Fetters vd. (1994), Tablolar 7–8](https://doi.org/10.1063/1.555949)

Örneğin IK1500-1’in [η] oranı 5306/667 ≈ 7,96 iken aynı fraksiyonun Rg küpleri oranı (506/228)³ ≈ 10,93’tür. İki ölçüm birbirini tamamlar. İntrinsik viskoziteden hidrodinamik genişlemeyi doğrudan geometrik yumak hacmine çevirmek ek bir model ve ortalama tanımları gerektirir.
