# Zincirin kuvvet alanı ve kısmi yükleri

Sonlu zincirin bütün kovalent bağları, uç grupları ve tanımlanan stereokimyası aynı molekül içinde ele alınır. Örnekte iki omurga bağlantısı olan azotlar tersiyer amid, H ile kapatılan uçtaki azot sekonder amiddir. Üç karbonil oksijeni su ile hidrojen bağı kabul edebilir; uçtaki N–H ayrıca donör olabilir. Başka bir uç grubu veya farklı zincir uzunluğu yeni bir yük dağılımı gerektirir.

1. İlk uygulamanın kaydettiği en düşük yakınsamış MMFF94s konformerinin tüm açık hidrojenleri ve koordinatları alınır. Kimyasal yapı, bağ sırası, atom sırası, koordinat birimi ve bildirilmiş stereokimya denetlenir.
2. AmberTools Antechamber atom tiplerini GAFF2 ile atar. SQM, bu tek başlangıç konformerinden gaz fazında yarıampirik AM1 geometri optimizasyonu ve Mulliken yüklerini hesaplar. Antechamber bağ yük düzeltmelerini uygulayarak AM1-BCC kısmi yüklerini verir. Bu işlem DFT içermez. Suda MD üretiminde sabit yükler kullanılır.
3. Çıktı yüklerinin toplamındaki küçük yazım veya yuvarlama sapması her atoma aynı miktar eklenerek giderilir. Bu düzeltmenin zincir toplamı 0,01 e veya atom başına 0,0005 e değerini aşarsa sonuç üretilmez. Pilotun işlenmemiş kısmi yükleri −0,001998 e toplam verdi. Elli atomun her birine +0,00003996 e eklenmesi MOL2 dosyasında toplamı sıfıra getirdi. Native Amber topolojisinde sayısal yazım sonrası toplam yaklaşık −9,0 × 10⁻⁹ e'dir.
4. Parmchk2 ek terimleri çıkarır. Eksik veya açıkça düzeltilmesi gereken bir parametre işlemi durdurur. Analog terim puanı 10 üzerindeyse bu sürüm otomatik MD'ye geçmez. Bu yazılım eşiğinin deneysel bir hata karşılığı yoktur. Pilot zincir için altı dihedral satır ve dört improper satırı eklendi. Üç improper analojisinin puanları 6, 3 ve 6 idi. Aktarılan terimler, amide ilişkin torsiyonlar dahil, tam biçimleriyle **solute.frcmod** dosyasındadır.
5. LEaP kovalent terimleri, Lennard–Jones terimlerini, dışlamaları ve 1–4 ölçeklerini native Amber topolojisine yazar. Yapı ve net yük OpenMM ile denetlenir. Pilotun native **solute.inpcrd** dosyası ilk uygulamanın koordinatlarını Antechamber çıktısının 0,001 Å hassasiyetinde taşır. SQM'nin yük hesaplamasındaki optimize koordinatları **sqm.out** içinde ayrıca bulunur.

Genel organik parametrelerin yapıyı kapsaması, PEtOx–su kuvvet alanının bu oligomer için deneysel doğrulaması değildir. Hidratasyon ve zincir boyutları AM1-BCC tek-konformer yük seçimine, amid torsiyonlarına, sabit yüklü su modeline ve örnekleme süresine duyarlıdır. Kaydedilen analoji puanları bu etkileri nicelleştirmez. Pilotun ölçülmüş bir öngörü hata payı yoktur.

## Yeniden üretim

Hesap paketi **pilot/parameters/parameterization.json** içinde kimyasal yapı, kaynak konformer kimliği, gerçek koordinat özeti, GAFF2 veri dosyaları ve program özetleri, analog parametre kayıtları, topoloji kontrolleri ile tüm çıktıların SHA-256 değerlerini taşır. Arşivde AmberTools çalıştırılabilir programları bulunmaz. Pilotun araçları Linux x86-64 için **acpype==2026.10.15** paketinden alınmıştır. GAFF2 veri dosyasının başlığında sürüm **2.2.30, Oct 2025**, SQM çıktısında **AMBER SQM VERSION 19** yazılıdır. AmberTools dağıtımının LICENSE dosyası, dat/leap kuvvet alanı parametrelerinin yazarlarınca kamu malı olarak yayımlandığını belirtir.

Aynı parametreler yeni bir sıcaklıkta kullanılabilir. Yeni kimyasal yapılar için RDKit, OpenMM ve gerçek bir AmberTools kurulumu gerekir. Uygulama önce **AMBERHOME**, ardından **PATH** ve ACPYPE'nin kullanılabilir Linux paketini arar. Uyumlu bir AmberTools kurulumuyla farklı platformlarda yeni parametre kaydı oluşturulabilir. Yazılım ve parametre özeti değişince yeni kimlik oluşturulur.

Proje dizininde aşağıdaki örnek, paketlenmiş verileri ilk uygulamanın kaydından yeniden hesaplar. **rebuild** seçeneği gerçek AmberTools programlarını tekrar çalıştırır. Tekrar üretimde prmtop başlığının oluşturma tarihi platforma ve çalıştırma saatine göre değişebilir.

```python
import json
from pathlib import Path
from parameters import parameterize

source = json.loads(Path("pilot/source_model.json").read_text())["energy"]
parameters = parameterize(
    source["canonical_smiles"], "cache/parameters-reproduction",
    source_record=source, rebuild=True,
)
print(parameters["parameter_id"], parameters["directory"])
```

Verilen kimyasal yapıda 2–100 ağır atom bulunmalı, olası stereo merkezleri ve çift bağ konfigürasyonu belirtilmelidir. İlk sürüm tek bileşenli, kapalı kabuklu, izotop etiketsiz, her atomu formal olarak yüksüz C/H/N/O yapıları işler. Parametreler veya model kontrolleri eksik kalırsa yeni çalıştırmanın sonuçları üretilmez.

## Kaynaklar

- Wang ve ark. (2004). General AMBER force field yönteminin temeli. [Journal of Computational Chemistry](https://doi.org/10.1002/jcc.20035)
- Jakalian, Jack ve Bayly (2002). AM1-BCC yük modelinin parametrelenmesi ve değerlendirilmesi. [Journal of Computational Chemistry](https://doi.org/10.1002/jcc.10128)
- Sousa da Silva ve Vranken (2012). ACPYPE yazılımı ve Antechamber araçlarının dağıtımı. [BMC Research Notes](https://doi.org/10.1186/1756-0500-5-367)
