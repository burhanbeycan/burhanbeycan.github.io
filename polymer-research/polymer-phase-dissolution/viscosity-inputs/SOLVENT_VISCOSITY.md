# Saf çözücünün dinamik viskozitesi

CoolProp 8.1.0 ve HEOS arka ucu, her saf sıvının yayımlanmış taşınım korelasyonunu çözelti numunesinin nominal sıcaklığında ve 101.325 Pa basınçta değerlendirir. Fazın sıvı olduğu yazılımla kontrol edilir. Hesaplanan dinamik viskozite Pa·s biriminden mPa·s birimine 1.000 ile çarpılarak çevrilir.

| Saf sıvı | T (K) | η₀ (mPa·s) | Viskozite korelasyonu |
|---|---:|---:|---|
| Benzen | 298,15 | 0,6021 | [Avgeri vd., 2014](https://doi.org/10.1063/1.4892935) |
| Toluen | 298,15 | 0,5522 | [Avgeri vd., 2015](https://doi.org/10.1063/1.4926955) |
| Toluen | 303,15 | 0,5205 | [Avgeri vd., 2015](https://doi.org/10.1063/1.4926955) |
| Siklohekzan | 307,65 | 0,7628 | [Tariq vd., 2014](https://doi.org/10.1063/1.4891103) |

Belirli koşulda çözeltinin mutlak viskozitesi η ≈ ηᵣη₀ olarak hesaplanabilir. ηᵣ kaynağının kapiler akış süreleri oranını kullandığı serilerde eşit yoğunluk yaklaşımı sonuca taşınır. Polimer kaydı başka bir basınç vermediğinde 101.325 Pa saf sıvı için kullanılan referans basınçtır. Copenhaver ölçümlerinin 30 ± 1 °C sıcaklık aralığı bulunduğundan referans değer 30 °C için nominaldir. Tablodaki basamak sayısı korelasyonun değerlendirme hassasiyetini yansıtır. Deneysel ve model belirsizliği ayrıca değerlendirilmelidir.

Hesabı tekrarlamak için Python ortamına `CoolProp==8.1.0` kurun ve `generate-solvent-viscosities.py` dosyasını çalıştırın. Oluşan `viscosity-solvent-data.js` her kayda sıcaklığı, basıncı, DOI'yi ve CoolProp uygulama bağlantısını ekler.
