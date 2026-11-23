# Polymer Research Studio 0.2.0

Burhan Beycan’ın polimer kimyası araştırmaları için hazırlanmış yedi ayrı kaynak proje. Uygulamalar tek tek açılabilir; aynı Python ortamına kurulduklarında ortak numune ve sonuç kayıtları üzerinden birlikte kullanılabilir. Arayüz İngilizcedir. Hesaplamalar bilgisayarınızda çalışır; API anahtarı gerekmez.

| GitHub için hazırlanmış proje | Bilimsel kapsam |
|---|---|
| `polymer-research-studio` | Akademik arayüz, ortak numune kaydı, sonuç geçmişi ve açıkça tanımlanmış veri aktarımı |
| `polymer-chemistry-lab` | Bağ bağlantıları belirli tekrar birimlerinden sonlu oligomerler, açık kopolimer dizileri, formül/kütle ve atom–bağ gösterimi |
| `polymer-synthesis-lab` | Fonksiyonel grup stokiyometrisi, Carothers, kondenzasyon kütle dengesi, ideal başlatılmış zincir hesabı ve anlık kopolimer bileşimi |
| `polymer-chain-lab` | Molar kütle dağılımları, uç grup düzeltmeli tekrar sayısı, Gaussian/WLC boyutları ve RIS konformasyon istatistiği |
| `polymer-transport-lab` | Tanımlı sınır koşullarında difüzyon, sabit geometrili sayısal doğrulama ve geçirgenlik |
| `polymer-network-lab` | Nötr ağ şişmesi, referans durumu açık termodinamik model ve denge geometrisi aktarımı |
| `polymer-electrolyte-lab` | Nernst–Einstein yaklaşımı, kolektif Onsager katsayıları ve referans çerçevesi tanımlı iyon taşınımı |

## Bilgisayarda kurulum

Python 3.10 veya üzeri gerekir. Yedi projeyi aynı klasöre çıkarın ve bu klasörde terminal açın. Önce sanal ortam oluşturun:

```bash
python -m venv .venv
```

Windows PowerShell: `.venv\Scripts\Activate.ps1`; macOS/Linux: `source .venv/bin/activate`. Ardından:

```bash
python -m pip install ./polymer-research-studio ./polymer-chemistry-lab ./polymer-synthesis-lab ./polymer-chain-lab ./polymer-transport-lab ./polymer-network-lab ./polymer-electrolyte-lab
polymer-studio
```

Ortak arayüz `http://127.0.0.1:8765` adresinde açılır. Bağımsız komutlar `polymer-chemistry-lab`, `polymer-synthesis-lab`, `polymer-chain-lab`, `polymer-transport-lab`, `polymer-network-lab`, `polymer-electrolyte-lab` şeklindedir. Tek laboratuvarı kullanmak için ana proje ile o laboratuvarı kurmanız yeterlidir. Python dağıtımının `polymer-research-workbench` adı, mevcut program arayüzleriyle uyumluluk için korunmuştur; GitHub ana proje adı `polymer-research-studio` olacaktır.

