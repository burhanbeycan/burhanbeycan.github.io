# Suda tek zincir için uzatılmış örnekleme

Bu protokol bir önceki 24 ps’lik üretimin son kaydedilmiş karesindeki **bütün atomları**, su moleküllerinin kimliklerini ve o andaki periyodik kutuyu devralır. İlk PEtOx–su arşivinde hızlar bulunmadığından önceki hızları gerektiren kesintisiz bir dinamik zaman serisi üretilemez. Her soy için Maxwell–Boltzmann hızları ayrı tohumla atanır, bağ kısıtları uygulanır, 12 ps NPT ile yeniden ısıl örnekleme yapılır. Uzatılmış bölümün ayrı üretim saati 0’dan başlar. Arşivlenen ilk kare 1 ps, son kare 120 ps’dedir.

## Yeniden çalıştırma

Polimer–çözücü uygulamasının açılmış klasöründen, uygulamanın Python bağımlılıklarını kurduktan sonra aşağıdaki komutları çalıştırın. İlk uygulamanın tekrar parametrelenmesi gerekmiyor.

```bash
for seed in 13031 27011 49033
do
    python hydration/experiments/extend_hydration.py \
        --base . --params pilot/parameters --source pilot/run --out hh-long \
        --seed "$seed" --production-ps 120 --discard-ps 12 --stride-ps 1
done
python trajectory_analysis.py --topology hh-long/topology.json --out hh-analysis \
    hh-long/seed_13031.npz hh-long/seed_27011.npz hh-long/seed_49033.npz
```

Windows’ta aynı Python komutu, `--seed` alanında her bir tohum belirtilerek üç defa çalıştırılabilir. Farklı üretim süresi veya örnekleme aralığı için yeni bir çıktı dizini verin. CROP uçlu zincirin ana koşusu önce kendi giriş yapısı ve kuvvet alanıyla hazırlanır; bu komutta ona ait parametre ve kaynak dizinleri kullanılır.

## Koşu kimliği ve örnekleme

- Kaynak NPZ, topoloji JSON, Amber `prmtop` ve `inpcrd`, kaynak kaydındaki SHA-256 değerleriyle karşılaştırılır. Yeniden oluşturulan sistemde tüm atom adları, elementleri, kütleleri ve bağ grafı doğrulanır. Kaynak Amber aktarım kodunun ve TIP3P parametrelerinin özetleri denetlenir.
- Üç soy, önceki üç ayrı başlangıçtan gelir. `SeedSequence([30971, seed])` yeni hız, Langevin integratörü ve Monte Carlo barostatı için ayrı tamsayı tohumları üretir. Yeni bir su kutusu yerleştirilmez. Başlangıçtaki koordinatlar kaydedilmiş `float32` değerlerinin çift duyarlığa dönüştürülmüş hâlidir. H içeren bağların kısıtları bir kez yeniden uygulanır.
- Sistem bir sonlu zincir ve üç noktalı rijit TIP3P sulardan oluşur. GAFF2/AM1-BCC kısmi yükleri ve parametreleri kaynak paketten gelir. Langevin-middle integratörü, 2 fs adım, 1 ps⁻¹ sürtünme, 298,15 K, 1 bar NPT, 25 adımda bir hacim denemesi, 0,8 nm gerçek uzay kesimi ve PME kullanılır. Bu dosyadaki değerler çalışmanın seçili protokolünü tanımlar. Herhangi bir bağımsız yeniden çalıştırmada gerçekleşen koşullar `seed_*.json` içinde tutulur.
- `seed_*.npz` içindeki `positions_nm`, `box_nm` ve `time_ps` bütün atomlar ve bağlı su hidrojenlerini içerir. Yoğunluk, sıcaklık, hücre potansiyel enerjisi ve zincirin en yakın periyodik kopyasına uzaklık aynı kareler için saklanır. Üretimdeki zaman aralığı, çok kısa su değişimi ve H-bağı kinetiği için bir zaman çözünürlüğü sınırı koyar.
- Çıktılar 12 karelik bölümlerle yazılır. Kesilmiş bir iş aynı dizinde sürdürülünce kaydedilmiş bölümler korunur. OpenMM bağlamının yerel kontrol noktası motor sürümüne ve sisteme bağlıdır. Yeniden başlatılan stokastik dizinin gelecekteki kareleri birebir aynı olmayabilir. `state_*.xml.gz` son kutu, koordinatlar, hızlar ve bağlam parametrelerini başka hesaplara aktarmaya yarar; ilgili Amber sistemi ayrıca gerekir.

Suyun seçilen sınırdan içeri girmesi, açık-su molekülünün belirli karbonil çevresinde bulunması ve su vericili hidrojen bağları bu fiziksel sistemin yerel gözlemlenebilir büyüklükleridir. Bu tek zincirli hücreden yığın polimerin çözünme hızı veya dengedeki çözünürlüğü doğrudan ölçülemez.
