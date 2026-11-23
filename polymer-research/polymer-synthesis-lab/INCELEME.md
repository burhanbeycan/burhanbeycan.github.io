# Reçeteden PEtOx molar kütlesi

Başlangıç beslemesi, sıcaklık ve süreden hesaplanan Mₙ, deneysel SEC sonucuna ne kadar yaklaşır?

## Uygulamada incele

İki PEtOx deney serisini ayrı çalıştır. Kesikli reaktör için dört zaman noktasını, akış reaktörü için dört beslemeyi incele. Özel CSV ile yayını, seriyi, derişimi veya beslemeyi dışarıda bırakabilirsin.

## Örnek sonuç

Kesikli seride sabit kinetik modelin medyan mutlak bağıl hatası %10,30 (P90 %16,75). Akış serisinde sırayla bir besleme dışarıda bırakıldığında sabit model %49,01, eğitim beslemelerine göre SEC ölçeği düzeltilmiş model %28,54 (P90 %46,66) verdi.

![İki PEtOx deney serisinde SEC molar kütlesine göre model hataları](docs/assets/benchmark-comparison.svg)

## Kimyasal yorum

SEC sonuçlarının PMMA kalibrasyonunu ve reaktör tipini kaydet. Gözlenen dönüşüm model girdisi değildir. Her kohort tek başlangıç monomer derişimi içeriyor. Yeni derişimde %10 altı hata hedefi açık.

[Bilimsel yöntem ve veri](docs/PETOX.md) · [Proje ve kullanım](README.md) · [Tek proje kaynak arşivi](../downloads/polymer-synthesis-lab-0.3.0.zip) · [Tüm projeler](../README.md)
