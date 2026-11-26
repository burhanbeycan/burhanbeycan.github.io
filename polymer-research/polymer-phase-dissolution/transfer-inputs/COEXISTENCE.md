# PS–siklohekzan eşfaz bileşimleri

## Birincil dış sınama

Nakata, Dobashi, Kuwahara, Kaneko ve Chu (1978), *Coexistence curve and diameter of polystyrene in cyclohexane*, Physical Review A 18, 2683. Numunenin ve yöntemlerin birincil bilimsel kaynağı [doi:10.1103/PhysRevA.18.2683](https://doi.org/10.1103/PhysRevA.18.2683).

Sayısal aktarım kaynağı Christian Wohlfarth’ın *CRC Handbook of Liquid-Liquid Equilibrium Data of Polymer Solutions* (2008) derlemesinde basılı sayfa 147, 1978NAK kaydı. Kayıttaki 0,0324 toplam polimer hacim kesri serisinin 20 satırı ayrı bir deneysel karşılaştırma kümesidir. “Sol phase” ile “gel phase” başlıkları sırasıyla polimerce fakir ve polimerce zengin sıvı faz olarak aktarılmıştır. Tabloda ikinci başlık ağ oluşumunu ölçmüyor. [Yayınevi](https://www.routledge.com/CRC-Handbook-of-Liquid-Liquid-Equilibrium-Data-of-Polymer-Solutions/Wohlfarth/p/book/9781032099675) · [Sayısal tablonun yer aldığı kopya](https://research-solution.com/uplode/books/book-37203.pdf)

Çözücü siklohekzan, polimer Pressure Chemical Company polistireni. Mₙ = 1.515.000 ve M𝓌 = 1.560.000 g/mol, buradan Đ = 1,02970297. Her satırda iki faz aynı sıcaklıktadır. Mutlak sıcaklık, tablodaki T꜀ = 303,65 K değerinden tablodaki sıcaklık farkı çıkarılarak hesaplanır. Kompozisyon, zengin faz ile fakir fazın her birindeki PS hacim kesridir. Hacme dönüşüm yapılmadı. Kaynak satırları mutlak sıcaklık kalibrasyonunun belirsizliğini veya ayrı dal hata çubuklarını vermiyor. T꜀ iki, sıcaklık farkı üç, hacim kesirleri dört ondalıkla yazılmıştır. Bunlar yazım hassasiyetidir. Basınç bu kayıtta belirtilmiyor.

Tablodaki seçili seri `tie_lines.csv` dosyasında bulunur. İkinci toplam bileşimden (0,0321) veri kullanılmamıştır. Bu seçim, tek malzemenin farklı başlangıç bileşimlerinin bağımsız dış numuneler gibi sayılmasını önler.

## Sayısallaştırmanın çapraz kontrolü

Hager, Anisimov, Sengers ve Gorodetskii (2002), *Scaling of demixing curves and crossover from critical to tricritical behavior in polymer solutions*, [arXiv:cond-mat/0206464v1](https://arxiv.org/pdf/cond-mat/0206464v1), Şekil 2 ve kaynaklar 42–44. Şekil, ikili PS–siklohekzan çözeltilerinin deneysel eşfaz bileşimlerini T (K) ve φ eksenlerinde yeniden çizer. Çizgi ve yıldızlar ayrıdır. Şeklin 1,56 milyon g/mol serisi Nakata 1978 ile aynıdır.

Dağıtımdaki `scripts/extract_hager_fig2.py` sabit sürüm PDF'nin SHA-256 değerini doğrular. Vektör eksenini sekiz yatay ve beş dikey ana işaretle kalibre eder, yalnız veri için çizilmiş dairesel işaret merkezlerini alır. Sıcaklık ve seri içinde aynen iki sembol, her daldan birer tane olan satırlar karşılaştırılır. Birden fazla olası çiftin aynı çizim satırına düştüğü durumlara eşleştirme yapılmaz. Bu işlemin 15 ayrı 1,56 milyon çiftini tablodaki karşılıklarıyla karşılaştırması en çok 0,00028 hacim kesri ve 0,023 K fark verir. Çizim kaynağında iki merkez arasında yaklaşık 0,238 PDF noktası sayısal adım vardır. Diğer seri için yaklaşık 0,0003 hacim kesri ve 0,03 K sayısallaştırma toleransı kullanılabilir. Bunlar deneysel hata çubuğu olarak kullanılamaz.

Çapraz kontrol betiği ayrıca M𝓌 = 200.000 g/mol örneği için 12 açıkça eşlenebilir çift çıkarabilir. Deney kaynağı Nakata, Kuwahara ve Kaneko, *J. Chem. Phys.* 62, 4278 (1975), [doi:10.1063/1.430349](https://doi.org/10.1063/1.430349). Birincil çalışmada iki eşfazın bileşimi diferansiyel refraktometreyle ölçülür. Derlemede etiketli aynı örneğin kesin Mₙ/Đ değerini doğrulayamadığımız için bu küme ek duyarlılık karşılaştırmasına ayrılmıştır. Önceki PP(P)DB serisindeki başka yayının Đ değeri buraya kopyalanamaz. 110.000 g/mol eğrisinin noktaları her iki fazda aynı sıcaklığa karşılık gelmediği için bu dosyalara alınmadı. Eşfaz verilerini okuma işlemi verileri açıklayan fizik modelinden bağımsızdır.
