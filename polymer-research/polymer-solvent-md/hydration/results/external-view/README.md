# Atomistik yörüngeleri VMD ve GROMACS’ta açma

Her dosya grubu bir çalışmanın **önceden seçilmiş en küçük numaralı tohumunun tamamını** içerir. Bütün su molekülleri ve suyun hidrojenleri saklanır. Kimlikler ardışık kareler arasında aynıdır.

| Dosya kökü | Başlangıç ve sistem |
|---|---|
| `hh-13031` | H/H uçlu tek PEtOx trimeri ve 500 su. İlk kare yeniden ısıl örnekleme sonrası 1 ps’de kaydedilmiştir. |
| `pair-73111` | Başlangıçta temas eden iki H/H uçlu PEtOx trimeri ve 700 su. 0 ps karesi tutucu potansiyel kaldırıldığında başlar. |

**VMD.** İstediğiniz kaydın `.pdb` dosyasını açın. Aynı moleküle `.xtc` dosyasını yörünge olarak ekleyin. PDB’nin ilk yapısı XTC’nin ilk karesiyle aynıdır. VMD’ye ayrı bir kare olarak yüklenmişse zaman analizinde bu yinelenen yapıyı çıkarın. Tek zincir için `resname POL`, çiftte `resname POA` ve `resname POB` seçilebilir. PDB’de sular `HOH` kalıntısıdır. Sadece yakın ortamın görüntülenmesi çözücünün kalanını simülasyon dosyalarından çıkarmaz.

**GROMACS.** `.gro` dosyası ilk tam sistemin atom sırasını, atom adlarını ve periyodik kutusunu taşır. `.ndx` içindeki gruplar standart bir tabanlı atom numaralarıyla oluşturulmuştur. `Water_O`, her suyun oksijenidir. `Amide_O`, `Amide_O_A` ve `Amide_O_B` karbonil oksijenlerini gösterir. Örneğin `gmx check -f pair-73111.xtc` kaydın kare ve zaman bilgisini okur.

Dosyalar birim hücrenin orijinal koordinatlarını içerir. Görüntüleyicide periyodik sınırı geçen su atlayabilir. Çift zıt yüzlerde görünürse fiziksel mesafe için en yakın periyodik kopya kullanılmalıdır. Uygulama tüm nicel yakınlık ve hidrojen bağı hesaplarında bu düzeltmeyi uygular.

GRO ve XTC koordinat hassasiyeti 0,001 nm’dir. Tüm XTC kareleri bağımsız olarak geri okunarak sıraları, zamanları, kutuları ve koordinatları bilimsel NPZ çıktılarıyla karşılaştırıldı. Her kaydın `.json` dosyası en büyük koordinat bileşeni farkını, atom ve su sayısını, zaman aralığını, giriş ve çıktı SHA-256 özetlerini verir. PDB/GRO konfigürasyonları bir MD kuvvet alanı içermez; yeniden hesap için ilgili Amber parametre paketi saklanmıştır.

### Bir başka koşuyu dışa aktarma

Uygulamanın kök klasöründen aşağıdaki örnek 27011 numaralı H/H koşusunu ayrı bir klasöre yazar.

```bash
python hydration/experiments/export_xtc.py --raw hydration/data/hh/run/seed_27011.npz --topology hydration/data/hh/run/topology.json --manifest hydration/data/hh/run/seed_27011.json --out my-view --name hh-27011
```
