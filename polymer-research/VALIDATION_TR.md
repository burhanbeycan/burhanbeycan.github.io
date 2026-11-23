# Sürüm 0.3.0 doğrulama kaydı

Yedi kaynak pakette 240 test çalıştırıldı (24 Kasım 2026). Ayrı oluşturulan yedi Python wheel paketi temiz bir hedef klasöre kuruldu. Kurulu uygulamaların yerel HTTP arayüzünden altı varsayılan hesap dâhil 49 örnek çalıştırıldı. Girdi, çıktı ve kaynak kayıtlarının dışa aktarım biçimi de kontrol edildi.

| Proje | Geçen test | Varsayılan dâhil kurulu paket örneği | Yöntem |
|---|---:|---:|---|
| [polymer-chain-lab](polymer-chain-lab/INCELEME.md) | 61 | 7 | [Yöntem](polymer-chain-lab/docs/PEtOx_VISCOMETRY.md) |
| [polymer-chemistry-lab](polymer-chemistry-lab/INCELEME.md) | 26 | 11 | [Yöntem](polymer-chemistry-lab/docs/SCIENCE.md) |
| [polymer-electrolyte-lab](polymer-electrolyte-lab/INCELEME.md) | 34 | 7 | [Yöntem](polymer-electrolyte-lab/docs/PEO_COMPARISON.md) |
| [polymer-network-lab](polymer-network-lab/INCELEME.md) | 31 | 7 | [Yöntem](polymer-network-lab/docs/WBPU.md) |
| [polymer-research-studio](polymer-research-studio/INCELEME.md) | 31 | ortak arayüz | [Yöntem](polymer-research-studio/docs/EXCHANGE.md) |
| [polymer-synthesis-lab](polymer-synthesis-lab/INCELEME.md) | 24 | 9 | [Yöntem](polymer-synthesis-lab/docs/PETOX.md) |
| [polymer-transport-lab](polymer-transport-lab/INCELEME.md) | 33 | 8 | [Yöntem](polymer-transport-lab/docs/PVA_WATER.md) |

**Arayüz kontrolü.** Gerçek API yanıtlarıyla JavaScript DOM etkileşimleri sınandı. PEtOx uçlarının çizilmesi, yakınlaştırma, aynı partiden iki ayrı numune türetme, aktif numuneyi değiştirme ve kaynağın doğru numunede görünmesi kontrol edildi. Masaüstü ve mobil tarayıcı yerleşimi görsel olarak doğrulanamadı.

**Deneysel doğrulama.** Viskozite ve şişme karşılaştırmalarında örnekler tek yayından geliyor. Sentezin derişimler arası Mₙ hedefi açık. PVA–su kalınlık karşılaştırması sentetik, nem karşılaştırması bir deneysel ara noktaya dayanıyor. Elektrolitin grafikten okunan katsayılarında okuma belirsizliği ayrıca veriliyor. Her projedeki yöntem sayfası kaynakları ve ölçüm koşullarını içerir.
