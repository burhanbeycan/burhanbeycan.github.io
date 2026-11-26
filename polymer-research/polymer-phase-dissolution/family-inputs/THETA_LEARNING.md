# θ koşullarında, önceden belirlenmiş bir referansla viskozite aktarımı

Hedef sonsuz seyreltmedeki intrinsik viskozitedir. Mol kütlesi g/mol, [η] mL/g cinsindedir. Her kaydın kimyası, çözücüsü, sıcaklığı, kütle tanımı ve birincil kaynağı korunur. Kaynak ölçümler ve sabit Ridge hesabı için aynı klasördeki `FAMILY_INTRINSIC_DATA.md` ve `THETA_MODEL.md` dosyalarına bakılabilir.

### Hedef aileden referans bulunmadığında

Önceki modelin dış katları korundu. Her turda bir omurga ailesinin tamamı saklandı. İki dış eğitim ailesindeki her kimya sırayla iç sınamaya ayrıldı. Sabit 0,5 kütle üssü ve yapı tanımlayıcısının cezasız eğimine uygulanan beş katsayı karşılaştırıldı (0, 0,25, 0,50, 0,75, 1). Aynı modelde 0,5 katsayısı önceki sabit Ridge cezasına λ=1 karşılık gelir. Her iç katın tanımlayıcı merkezi ve eğimi, o katın eğitim gözlemleriyle yeniden hesaplandı. İç sınamada aile, kimya, çalışma ve numune için önceki ağırlıklarla log-hata karesi ortalandı. En küçük hatadaki katsayı seçildi. Eşitlikte daha küçük eğim kullanıldı.

Üç aileye eşit ağırlık verildiğinde eski ve uyarlanmış modeller için 56 hedefte MAPE sırasıyla %12,71 ve %13,00 oldu. İki eğitim ailesinin az sayıdaki kimyası, düzenlileştirmeyi öğrenmek için sınırlı çeşitlilik sağlıyor. Bu sonuç sabit Ridge'in kullanımını destekliyor.

### Bir literatür fraksiyonu varsa

Referans aynı kimya, çözücü, bildirilen nominal θ sıcaklığı ve aynı mol kütlesi tanımına sahip ölçümler arasından seçilir. 100.000 g/mol'e log ölçekte en yakın fraksiyon kullanılır. Eşitlikte numune kimliğinin sırası belirleyicidir. Seçim [η] değerini kullanmaz. Polidienlerin derlemedeki Mabs değerleri, Mw olarak etiketlenmez. PI'nin yaklaşık sıcaklığı ve mikroyapıdaki olası dağılım, kaynak dosyasında işaretlidir.

Üç karşılaştırmanın ilki hedef aileyi görmemiş eski yapı modelidir. İkinci model tek bir referansı θ güç yasasıyla aktarır.

    η_phys(M) = η_ref × √(M/M_ref)

Üçüncü model dış eğitimdeki kimyalar sırayla saklandığında gözlenen log-viskozite hatalarından ne kadar düzeltme uygulanacağını öğrenir. Her eğitim kimyası için kimya-dışı Ridge'in log-hataları hesaplanır. Kimya başına ağırlıklandırılmış ortalama hatanın karesinin ortalaması B, aynı kimya içindeki hata varyanslarının ortalaması W olarak hesaplanır. Aile ve kimya ağırlıkları korunur. Varyansın paydası etkili örneklem sayısına göre düzeltilir. w=B/(B+W), sıfır payda için w=0 alınır. B saf bir popülasyon varyansı değildir. Bu, tek ölçümlü aktarımı dengeleyen moment tabanlı bir kestirimdir; tahmin aralığı oluşturmaz.

    log η_adapt(M) = log η_ridge(M) + w × [log η_ref − log η_ridge(M_ref)]

Beş referans gözlemi 56 kayıttan ayrıldı. Kalan 51 hedefte aynı gözlemler karşılaştırıldığında aile ortalaması MAPE eski Ridge için %12,36, yalın fiziksel referans için %2,78, öğrenilmiş kısmi düzeltme için %3,10'dur. Referansların çoğu hedefleriyle aynı çalışmadan gelir. Bu karşılaştırmadaki kazanım, hedef kimyaya ait bir deneysel bilginin eklenmesini de içerir. Saf moleküler yapıdan elde edilen tahmin becerisi, sadece ilk 56 örneklik referanssız protokolde sınanır.

### Referans yayınının dışına aktarım

İki farklı çalışmanın PS–siklohekzan (34,5 °C) verilerinde her çalışmanın tüm satırları sırayla hedef grubudur. Referans yalnızca diğer çalışmanın aynı türdeki satırlarından seçilir. Veri tabanında bu bağımsız çalışma karşılaştırması için yeterli kaydı bulunan tek kimya PS'dir. Çalışmalara eşit ağırlıkla MAPE, eski modelde %12,34, yalın referansta %5,62 ve kısmi düzeltmede %5,58'dir. Buradaki %0,04 puan fark iki çalışma ile güçlü bir AI katkısı göstermek için yeterli değildir. Çalışmaların farklı mol kütlelerini örneklemesi, kaynak ve kütle aralığının etkisini bir arada tutar.

Arayüzdeki bilinen bir hedef numune kendi kalibrasyon referansı olarak seçilmez. Seçilen numunenin kendisi 100 kg/mol'e en yakın fraksiyonsa, arayüz o kimya için bir sonraki uygun fraksiyonu kullanabilir. Araştırma sınamasında ilk beş standart performans paydasından tamamen çıkarılır. `theta-ps-frozen.json`, yeni bir PS yayınının değerlendirilmesi için dış hedefler incelenmeden dondurulan önceki model, katsayılar ve referans kimliğini içerir.

### Yeniden hesaplama

```sh
node scripts/build-theta-learning.cjs dist dist/family-inputs
node tests/test-theta-learning.cjs dist
```

Düz çevrimdışı pakette `dist` yerine `.` kullanılabilir. `theta-learning-results.json` kaynak satırlarının tahminlerini, katsayıları ve çalışmaya/aileye göre hataları içerir. Bağıl hata yüzdeleri mutlak tahmin hatasının ölçüme bölünmesiyle hesaplanır.
