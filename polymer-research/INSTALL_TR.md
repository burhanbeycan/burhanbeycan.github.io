# Polymer Research Studio 0.4.0 kurulumu

**Gereksinim:** Python 3.10 veya üstü. Yedi projeyi bir arada içeren [ZIP arşivini](downloads/polymer-research-suite-0.4.0.zip) aç ve terminalde `polymer-research` klasörüne gir. Git kullananlar aynı klasöre şu şekilde ulaşabilir.

```bash
git clone --single-branch --branch polymer-research-review-v0.2.0 https://github.com/burhanbeycan/burhanbeycan.github.io.git
cd burhanbeycan.github.io/polymer-research
```

Sanal ortam oluştur.

```bash
python -m venv .venv
```

macOS ve Linux için `source .venv/bin/activate`, Windows PowerShell için `.\.venv\Scripts\Activate.ps1` çalıştır. Sonra yedi klasörden paketleri kur.

```bash
python -m pip install ./polymer-research-studio ./polymer-chemistry-lab ./polymer-synthesis-lab ./polymer-chain-lab ./polymer-network-lab ./polymer-transport-lab ./polymer-electrolyte-lab
polymer-studio
```

Tarayıcıda `http://127.0.0.1:8765` açılır. Sol menüden laboratuvarı ve **Reference scenario** menüsünden örneği seçip **Run calculation** düğmesine bas. **Model & references** varsayımları ve kaynakları, **Material record** fiziksel numuneleri, **Workflow** hesap geçmişini gösterir. Oturumu saklamak için **Save workspace**, tekrar açmak için **Import record** kullan.

Tek bir uygulama için Studio klasörünü ve istediğin laboratuvarı aynı ortama kur. Örneğin `python -m pip install ./polymer-research-studio ./polymer-chain-lab` ardından `polymer-chain-lab`. İlk kurulum için pip bağımlılıkları indirir. Kurulumdan sonra hesaplar yerel olarak çalışır.

## Yeni pilotu inceleme

Studio’da **Research pilot** sekmesi çalıştırılmış 16 kısa yörüngenin sonuçlarını kurulu paketten okur. **Starting pose**, **Partner sequence** ve **Acceptor site** ile koşulları değiştir. Sıfırdan çalıştırmak için [LDI pilot dosyasını](polymer-research-studio/examples/ldi-blend/README.md), PEtOx–su kümesi için [Chain Lab yöntemini](polymer-chain-lab/docs/PEtOx_MD.md) aç. Yeni film simülasyonlarına başlamadan önce numuneye ve kuvvet alanına uygun parametreler ile eşleştirilmiş deneysel ölçüm gerekir.
