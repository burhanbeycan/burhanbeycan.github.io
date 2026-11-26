# İyi çözücüde seyreltilmiş zincirlerin intrinsik viskozitesi

## Bir yayımdan diğerine parametre aktarımı

Çalışma lineer, nötr ve yüksek 1,4 içerikli polidienleri ve poliizobütileni nominal 25 °C’de siklohekzanda inceler. Hedef seyreltik limitte [η] (mL/g), giriş mol kütlesi (g/mol) ve kayıtlı kimyadır. Uygulama, farklı numunelere ait uygun referansları deney kataloğundan seçer.

Kimya c dışarıda bırakılarak diğer kimyalara ln[η] = b + a ln M modeli uydurulur. Her eğitim kimyası eşit ağırlık alır. Tek referanslı yöntemde eğim korunur ve referans numunenin ln[η] değeriyle kesişim yeniden kestirilir.

    [η]tahmin = [η]referans × (Mhedef / Mreferans)^a

Bu işlem bir parametreyi az veriyle uyarlayan log-doğrusal transfer modelidir. Hedef kimyanın yayımlanmış bir [η] ölçümü bu yolun açık ek girdisidir. Kalibrasyon, yüksek derişim veya kayma mekanizması hakkında bilgi sağlamaz. Hedef numunenin Rg değerine gerek yoktur.

**Referans seçimi.** Ortam, nominal sıcaklık, kimya ve katalogda raporlanan M türü eşleştirilir. Hedef çalışmanın DOI’sine ait tüm kayıtlar çıkarılır, kalanlar içinde |ln(Mhedef/Mreferans)| asgari olan seçilir. Eşit mesafede seçim kararlı numune kimliğine göre yapılır. Kimyasal mikroyapının, dağılım genişliğinin ve gerçek sıcaklığın kayıtlardaki ayrıntı düzeyiyle sınırlı eşleşme söz konusudur. Derleme kimi benzer iyi çözücü ölçümlerini küçük sıcaklık farklarıyla birlikte değerlendirir. PBD ve PI için M_table, kütlenin Fetters derlemesindeki M tanımına göre alındığını gösterir. Bu etiket tüm fraksiyonlar için aynı deneysel ortalama veya mol kütlesi yönteminin kanıtı değildir. PIB serileri M𝑤 kullanır.

**Yeni kullanıcı numunesi.** Kimya ve M ile tahmin yapılır. M aynı kimyadaki bir yayımlanmış numunenin değerine tam eşleşiyorsa o değeri bildiren çalışmalardaki tüm satırlar referans havuzundan çıkarılır. Kullanıcı yayındaki bir numuneyi seçerse gerçek DOI ile ayrım yapılır. Kayıt kalmadığında uygulama tek referanslı sayı döndürmez. Her numune için kullanılan referans, sayısal M aralığı ve kalibrasyon üsleri izlenebilir. Bilinmeyen bir numunenin yalnız M değerinden üretim tarihçesi veya aynı partiden gelme durumu kesinleştirilemez.

## Eşleşen 18 hedef ve iki kimya

22 özgün M–Rg–[η] kaydından iki ayrı makalede ölçüm bulunan PI ve PIB kimyaları 18 hedef sağlar. Dört PBD kaydı aynı yayındadır. Bu dört kayda başka PBD kaydı atayıp aynı yayını sınamanın iki tarafına koymak yerine, bunlar bu aşamada sınama dışında bırakıldı. Her hedef yayını tek seferde ayrılmıştır. Geriye kalan aynı kimyadaki yayından seçilen en yakın referans, özellikle tek kayıt olduğunda birden fazla hedef için kullanılabilir. Dört hedef yayını ve iki kimya vardır.

| Hedef kimya | n | Diğer kimyaların mol kütlesi modeli | Bir referansla aynı üssün aktarımı | Rg³/M modeli | Rg³/M için Bayesçi düzeltme |
|---|---:|---:|---:|---:|---:|
| PI | 6 | %8,68 | %8,38 | %49,57 | %40,73 |
| PIB | 12 | %27,85 | %6,43 | %11,71 | %18,37 |
| Kimyalara eşit ağırlık | 18 | %18,26 | %7,41 | %30,64 | %29,55 |

Her sayı aynı satırlardaki ortalama mutlak yüzde hatadır. Ölçümler ve kaynaklar önceki sürümden biliniyordu. Karşılaştırma retrospektiftir. Ekranda en düşük hatalı sürüm gösterilirken öteki iki hipotezin sonuçları da korunur. İki kimya yeni ailelere ait hassas bir genelleme hata aralığı belirlemeye yetmez.

## Yeni hedef çalışma

Colby, Fetters ve Graessley’nin makalesinde seçili dar dağılımlı, düşük vinil içerikli polibütadienlerin 25 °C'deki çözelti viskozimetresi ayrıca raporlanmıştır. Tablo II'nin siklohekzan sütunundaki sekiz dolu hücre alındı. Özgün [η] (dL/g) değerleri 100 ile çarpılarak mL/g elde edildi. Polimerlerin eriyik viskozitesi bu çalışmanın hedefi değildir. Tablo II'nin M=70.900 ve [η]=1,424 dL/g basamakları gibi değerlerde birincil kaynaktaki hassasiyet kullanıldı. Mol kütlelerinin tayininde ışık saçılması ve SEC, ortalamanın hesaplanmasında numune ve kütle aralığına bağlı olarak farklılık gösterir. Mikroyapının endüstriyel PBD çeşitlerine aktarılabilirliği bu sonuçtan bilinemez.

Bu hedef çalışmanın hiçbir verisi yöntem seçimine, üs uyumuna veya referans havuzuna girmedi. Kütle modelinin üssü 18 PI ve PIB kaydından hesaplandı. Referans için Roovers–Martin makalesindeki derleme tarafından **lineer PBD** grubuna alınmış dört kayıt kullanılabildi. En yakın log mol kütlesi ölçütü sekiz hedefte bu kayıtlardan iki farklısını seçti. Referans havuzu 114–361 kg/mol, dış hedefler 20–3.740 kg/mol aralığındadır. Sekiz hedefin altısı bu referans aralığının dışındadır. Aynı sekiz ölçümde MAPE %27,57 → %4,61, MAE 89,51 → 18,12 mL/g oldu. Bu sonuç bilinen kimyada başka bir çalışmaya aktarımı sınar. Yeni kimya veya yeni çözücüye aktarım için ayrı değerlendirme gerekir. Tüm numuneler tek bir çalışmadan gelir. [Colby ve ark., 1987, Tablo II ve deney yöntemi](https://doi.org/10.1021/ma00175a030) · [Roovers ve Martin, 1989](https://doi.org/10.1002/polb.1989.090271209) · [Fetters ve ark., derleme, Tablo 2](https://srd.nist.gov/jpcrdreprint/1.555949.pdf)

## Model birleştirme ve Bayesçi düzeltme

Birleşim için log tahmin = (1 − w) ln [η]kütle + w ln [η]boyut. İç katlarda eğitim kimyalarının tamamı birer kez tutulur. Kısıtlı en küçük kareler katsayısı 0 ≤ w ≤ 1 aralığındadır. Her kimyanın toplam ağırlığı aynıdır. Dış sınamanın [η] değerleri iç uyuma girmez. Yalnız tek kimya eğitimde kaldığında iç uyum mümkün değildir. Hedef çözücü veya nominal sıcaklık değiştiğinde ve iç kimya ayrımı mümkün olmadığında eski boyut modeli aynen korunur. Önceki hidrodinamik başarısızlıklar bilinirken tasarım yapılmıştır. İçte model öğrenilmesi bu tarihsel karşılaştırmayı ileriye dönük deneye dönüştürmez. [Breiman, 1996](https://doi.org/10.1023/A:1018046112532)

Üç kimyanın her biri için tüm kayıtlar dışarıda tutulduğunda aile dengeli MAPE yalnız mol kütlesinde %20,79, boyutta %29,18 ve birleşimde %34,51 olur. Boyut ağırlığı PB’de 0,386, PI’da 1 ve PIB’de 0 olarak bulunur. Aynı solventte yalnız iki eğitim kimyasıyla iç model seçiminin duyarlılığı görülür. Dış çözücüde çözeltinin ortak kütle yasası doğrulanmadığından, birleşim eski boyut modelinden ayrılmaz.

Bayesçi seçenek log hidrodinamik oran için diğer kimyalar arası varyans ile kimya içi varyansın parametrelerini yöntemin momentleriyle öğrenir. Her kimyada örnek varyansı s² ve örnek sayısı n ile kimya ortalamaları arasındaki örnek varyansından ortalama(s²/n) çıkarılarak τ² ≥ 0 bulunur. Tek bir referans ölçümü için ampirik Bayes ağırlığı w = τ² / (τ² + ortalama(s²)) olur. ln U = (1 − w) ln Ueğitim + w ln Ureferans, U = [η] M / ((4π/3) NₐRg³). İç varyansın kütleye göre sabit olduğu ve log ölçekte normal yaklaşık dağılım varsayımları vardır. İki eğitim kimyası, kararlı bir τ² kestirimi için sınırlıdır. PI ve PIB sonuçları ters yönlüdür. Boyut yolu için hedef Rg, tercihen bağımsız saçılma ölçümünden gelmelidir. [Polimerde viskozite yarıçapı, çözücü kalitesi ve U](https://doi.org/10.1021/ma500960f)

## Yeniden üretim

```sh
node scripts/build-hydro-stack.cjs dist
node scripts/hydro-protocols.cjs dist/hydro-stack-results.js hydro-protocols.json
node tests/test-hydro-improvement.cjs
```

Düz dağıtımda ilk iki komutta `dist` yerine `.` kullanılır. Testler tüm 26 referanslı tahminin ekrandaki hesapla eşleşmesini, dış çalışmanın sonuç etiketleri değiştiğinde tahminlerin sabit kalmasını, kimya ve DOI izolasyonunu ve negatif bulguları denetler.
