# Saf sıvı referansı · Su

Bu bağımsız modül periyodik bir kutuda 216 rijit TIP3P su molekülünün klasik NPT moleküler dinamiğini çalıştırır. Hedef koşul 298,15 K ve 1 bardır. Başlangıç geometrisini ve TIP3P kuvvet alanını OpenMM sağlar. Potansiyel, ayrıca üretilen izole MMFF94s çözücü kaydından bağımsızdır.

## Çalıştırma

Python ortamında `python -m pip install -r requirements.txt` komutuyla bu modülün bağımlılığını yükleyin. `python run_water.py all --out yeni_kosu` komutu üç replikayı ve ortak özeti üretir. Replikaları ayrı süreçlerde yürütmek için `python run_water.py run --seed 13031 --out yeni_kosu` komutunu 13031, 27011 ve 49033 tohumları için çalıştırın. Ardından `python run_water.py summarize --out yeni_kosu` komutunu verin.

Tüm replikalar aynı başlangıç kutusundan başlar. Başlangıç hızları, Langevin dinamiği ve Monte Carlo barostatı ayrı tohumlar kullanır. Her replika 12 ps başlangıç adaptasyonundan sonra 32 ps üretim içerir. Örnekleme aralığı 0,2 ps, zaman adımı 2 fs'dir. Uzun menzilli elektrostatik için PME, gerçek uzay ve Lennard–Jones etkileşimleri için 0,8 nm kesim kullanılır. Kuvvet alanı dosyasının özeti, kutu içeriği, hesap motoru sürümü ve diğer koşullar sonuç kaydındadır.

## Dosyaları kullanma

- `water_npt.csv` her replika ve aşama için hacmi, yoğunluğu, kinetik sıcaklığı ve kutunun potansiyel enerjisini içerir. Yoğunluk kutudaki toplam kütlenin ölçülen periyodik kutu hacmine bölünmesiyle bulunur.
- `water_npt.json` ortalamaları, tohumları, zaman serilerini ve dosya özetlerini taşır. Yoğunluk özetlerinde yalnız üretim bölümü kullanılır. Replika ortalamalarının örnek standart sapması ve her replikanın dört adet 8 ps blok ortalaması ayrı raporlanır. Aynı üretim karelerinin hücre potansiyel enerjisi 216 su molekülüne bölünür. Bir mol H₂O başına kJ olarak verilen ortalama bu kuvvet alanının referansına ve bu protokole özgüdür.
- `water_*_state.xml.gz` her replikanın son periyodik kutusunu, tam duyarlıklı atom koordinatlarını ve hızlarını içerir. Uyumlu TIP3P topolojisine sahip yeni bir koşu için başlangıç durumu olarak okunabilir. Modeller ve çözünen bileşen değiştiğinde yeni sisteme uygun hazırlama, dengeleme ve örnekleme gerekir. Rastgele sayı üretici durumunu da koruyan birebir simülasyon devamı için OpenMM checkpoint dosyası gerekir.

Bu kısa, sonlu-kutu pilotunda tek bir su modelinin yoğunluğu örneklenir. Daha büyük kutu ve daha uzun süre karşılaştırması, dengeleme penceresine duyarlılık ve uyumlu deneysel referansla doğrulama ayrı çalışmalardır. Polimer–su etkileşim enerjisi, difüzyon, çözücü viskozitesi veya solvatasyon serbest enerjisi bu kayıttan çıkarılmaz.

Model ve yöntem belgeleri: [OpenMM ForceField](https://docs.openmm.org/latest/userguide/application/02_running_sims.html), [su ile periyodik sistem kurma](https://docs.openmm.org/latest/userguide/application/03_model_building_editing.html).

## GROMACS koordinat dosyası

`water_13031.gro` bir koşunun son kutusunu, aynı atom sırasıyla ve periyodik kutu boyutlarıyla taşır. Diğer koşu için `python export_gro.py --seed 27011 --out water_27011.gro` komutu kullanılabilir. Dosyada uzunluk nm, hız nm/ps cinsindendir. GROMACS ile üretim hesabı için aynı TIP3P tanımı ve uyumlu polimer parametreleriyle topoloji hazırlanır. Kaynak MD koşuları OpenMM ile yürütülmüştür.
