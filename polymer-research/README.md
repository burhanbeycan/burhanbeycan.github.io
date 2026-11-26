# Polymer Research Studio 0.4.0

LDI–PEG/PCL esaslı üretan zincirlerini, PEtOx hidrolizini ve diziye bağlı hidratasyonu yedi ayrı projede inceleyebilirsin. İlk olarak Studio’daki **Research pilot** sekmesinde başlangıç pozunu değiştir. İki nötr model dimer ile etil-LDI diüretan arasındaki kısa, çözücüsüz klasik MD koşularını başlangıç bazında göreceksin.

| Proje | Çalışan örnek | Sonucu yorumlarken |
|---|---|---|
| [Studio ve WBPU/PEtOx pilotu](polymer-research-studio/INCELEME.md) | 2 dimer × 2 başlangıç kuralı × 4 tohum; N–H···O temaslarını karşılaştır | Yerel molekül modeli; kısa yörüngeler başlangıç konumuna duyarlı |
| [Network Lab](polymer-network-lab/INCELEME.md) | LDI eşdeğerleri, BDO/gliserol, periyodik film temasları ve kuru film FTIR karşılaştırması | Gerçek reçeteyle sayısal örnek ayrılmıştır; eşleşmiş film FTIR ve MD karşılaştırması için kullanıcı verisi gerekir |
| [Chemistry Lab](polymer-chemistry-lab/INCELEME.md) | 2026 NMR serisiyle PVA dizi blokluluğu; EtOx/EI dönüşümü, protonlanma ve karşı iyon | Altı PVA numunesi tek kaynaktan; EI mikrohalini kullanıcı tanımlar |
| [Chain Lab](polymer-chain-lab/INCELEME.md) | 12 su moleküllü PEtOx kümesinde iki zincir uzunluğu ve iki uç için 12 gerçek hesap | Kontrol edilen nicelik yerel su komşuluğudur |
| [Synthesis Lab](polymer-synthesis-lab/INCELEME.md) | Dört akış ölçümünden altı aday arasından yeni derişim ve besleme seçimi | Önerilen aday henüz ölçülmemiştir |
| [Transport Lab](polymer-transport-lab/INCELEME.md) | Aynı numunede bileşim, NMR ve hidratasyon bilgisiyle üç dış sınama | Üç kaynağın taramasında tam eşleşmiş deney serisi bulunamadı |
| [Electrolyte Lab](polymer-electrolyte-lab/INCELEME.md) | Li–eter ve Li–TFSI koordinasyonu, temas süreleri ve kolektif iletkenlik | Yayımlanmış MD karşılaştırması ile kullanıcının yüklediği yörüngeler ayrı değerlendirilir |

[Kurulum](INSTALL_TR.md) · [Yedi projenin kaynak arşivi](downloads/polymer-research-suite-0.4.0.zip) · [Doğrulama](VALIDATION_TR.md)

## Yapı ve konformer enerjisi

[Kendi polimerini oluştur ve konformer enerjilerini incele](polymer-structure-energy/README.md). Bu bağımsız uygulama, tanımlanan sonlu zincirin kimyasal kaydını ve izole konformer enerjilerini diğer modüllerin yeniden kullanabileceği biçimde saklar. Çözücü kayıtlarında yöntem ve hesap kapsamı izlenebilir.

[Bağımsız uygulamayı indir](downloads/polymer-structure-energy-1.0.0.zip).

## Polimer–su moleküler dinamiği

[PEtOx’un suda zincir boyutunu ve karbonil çevresindeki su dağılımını incele](polymer-solvent-md/README.md). Bu uygulama sonlu zincirin kaydını yeniden kullanır. Üç bağımsız kısa yörüngeden jirasyon yarıçapını, karbonil oksijeni ile su oksijeni arasındaki radyal dağılımı ve seçilen uzaklık içindeki ortalama su sayısını hesaplar.

[Grafikler](polymer-solvent-md/pilot/overview.png) · [Çevrimdışı etkileşimli raporu indir](https://github.com/burhanbeycan/burhanbeycan.github.io/raw/refs/heads/polymer-research-review-v0.2.0/polymer-research/polymer-solvent-md/pilot-review.html?download=1) · [Uygulama, çevrimdışı rapor ve ham yörüngeler (ZIP)](downloads/polymer-solvent-md-1.0.0.zip)

### PEtOx zincirinin suyla hidratasyonu

Modelin 1.1 sürümü, H/H uçlu PEtOx'un örneklemesini uzatır, kimyasal uç gruplarını karşılaştırır ve birbiriyle temas eden iki oligomerin sudaki hareketini gerçek atomistik karelerden izletir. Koşu girdileri, ham koordinatlar ve analizler birlikte sunulur.

[Hidratasyon çalışması ve sonuçları](polymer-solvent-md/hydration/README.md) · [Uygulamayı ve ham verileri indir](downloads/polymer-solvent-md-1.1.0/README.md)

## Polimer–çözücü faz davranışı ve çözünme

[Yapı, faz dengesi ve kinetik uygulaması](polymer-phase-dissolution/README.md) · [İki PEtOx zinciri arasına giren suyu izle (MP4)](polymer-phase-dissolution/PEtOx_Iki_Temas.mp4)
