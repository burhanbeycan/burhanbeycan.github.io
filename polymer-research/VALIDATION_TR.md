# Doğrulama kaydı · 0.4.0

Yedi kaynak pakette 298 yerel test geçti. Yedi ayrı wheel paketini kurduktan sonra altı varsayılan hesap dâhil 63 örnek, yerel HTTP uçları üzerinden çalıştırıldı ve sonuç belgeleri doğrulandı.

| Proje | Geçen test | Kurulu paket örnekleri |
|---|---:|---:|
| [polymer-chain-lab](polymer-chain-lab/INCELEME.md) | 68 | 8 |
| [polymer-chemistry-lab](polymer-chemistry-lab/INCELEME.md) | 34 | 13 |
| [polymer-electrolyte-lab](polymer-electrolyte-lab/INCELEME.md) | 46 | 10 |
| [polymer-network-lab](polymer-network-lab/INCELEME.md) | 42 | 10 |
| [polymer-research-studio](polymer-research-studio/INCELEME.md) | 34 | ortak arayüz |
| [polymer-synthesis-lab](polymer-synthesis-lab/INCELEME.md) | 31 | 11 |
| [polymer-transport-lab](polymer-transport-lab/INCELEME.md) | 43 | 11 |

## Uygulama ve model kontrolleri

Arayüz olayları gerçek API yanıtlarıyla bir DOM ortamında sınandı. PEtOx uçlarının çizilmesi ve yakınlaştırılması, fiziksel numune ayrımı, kaynak kayıtları, LDI pilotunda başlangıç konumu ve dimer değişimi, ayrıca eksik simülasyon koşullarının işaretlenmesi çalıştı. Tam masaüstü ve mobil tarayıcı görünümü bu kontrole dahil değildir.

LDI için 16 çözücüsüz ve PEtOx için 12 adet 12-su-moleküllü kısa MMFF94 yörüngesi çalıştırıldı. Bunlar hesaplanmış sonuçlardır. Kütlesel denge, FTIR bant alanı, gerçek film afinitesi ve elektro-eğrilmiş lif performansı doğrulaması sağlamazlar. PVA NMR ve PEO koordinasyon karşılaştırmaları kaynakların yayımlanmış değerlerini kullanır. Aktif öğrenme seçimleri ve sentetik denetimler her modülde işaretlenmiştir. Bu ortamda periyodik GROMACS üretim koşusu gerçekleştirilmedi.
