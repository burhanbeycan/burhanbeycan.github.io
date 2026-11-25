"""Build a portable Turkish interpretation from the frozen computed results."""
import argparse
import base64
import csv
from html import escape
import json
from pathlib import Path


def tr(x, digits=2):
    if x is None:
        return "—"
    return f"{x:.{digits}f}".replace(".",",")


def read(path):
    return json.loads(path.read_bytes())


def table(headings, rows):
    heads="".join(f"<th>{escape(str(x))}</th>" for x in headings)
    body="".join("<tr>"+"".join(f"<td>{escape(str(x))}</td>" for x in row)+"</tr>" for row in rows)
    return f"<div class=scroll><table><thead><tr>{heads}</tr></thead><tbody>{body}</tbody></table></div>"


def picture(path, alt, style=""):
    uri=base64.b64encode(path.read_bytes()).decode()
    return f'<figure class="{style}"><img src="data:image/png;base64,{uri}" alt="{escape(alt)}"></figure>'


def section(name, topic, title, inner):
    return f'<section class="sheet" id="{name}"><p class="topic">{topic}</p><h2>{title}</h2>{inner}</section>'


STYLE="""
:root{color-scheme:light;--ink:#132d38;--muted:#435965;--teal:#006f68;--line:#c9d7d7;--paper:#fff;--wash:#edf4f3}
*{box-sizing:border-box}body{margin:0;background:#eaf0ef;color:var(--ink);font:16px/1.48 system-ui,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif}
header{color:white;background:#14323d;border-bottom:4px solid #4da391}header>div,nav>div{max-width:980px;margin:auto;padding:12px 26px}header span{font-family:Georgia,serif;font-size:21px}nav{background:#fff;border-bottom:1px solid var(--line)}nav>div{display:flex;flex-wrap:wrap;gap:8px 20px}
nav a,a{color:var(--teal);font-weight:650;text-underline-offset:.19em}
main{max-width:980px;margin:25px auto 50px;padding:0 18px}
.sheet{margin:20px 0;background:white;border:1px solid var(--line);border-radius:12px;padding:clamp(20px,4.4vw,50px);scroll-margin-top:15px}
.topic{color:var(--teal);margin:0 0 4px;font-size:.92rem;font-weight:700}h2{font:clamp(26px,4.2vw,42px)/1.11 Georgia,serif;letter-spacing:-.02em;margin:0 0 14px}h3{font-size:1.1rem;margin:1.2em 0 .4em}p{margin:.65em 0}figure{margin:14px 0}img{max-width:100%;display:block;height:auto;margin:auto}.scroll{overflow:auto}table{border-collapse:collapse;font-size:.94rem;width:100%;margin:10px 0}th,td{padding:8px 10px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}th{background:#e5efed;color:#103d43}td:first-child{font-weight:600}ul{padding-inline-start:1.35rem}li+li{margin-top:.35rem}.note{background:var(--wash);padding:11px 14px;border-left:3px solid #4a8e80}.sub{color:var(--muted);font-size:.92rem}.small{font-size:.88rem}.button{display:inline-block;border-radius:5px;text-decoration:none;color:white;background:var(--teal);padding:8px 14px}
@media(max-width:600px){main{padding:0 9px;margin-top:12px}.sheet{border-radius:8px;margin:10px 0;padding:19px 14px}th,td{padding:6px 7px}.scroll table{font-size:.85rem}header>div,nav>div{padding:10px 15px}}
@page{size:A4;margin:13mm 14mm 15mm}
@media print{body{background:white;font-size:9.8pt;line-height:1.3}header,nav,.button,.screen{display:none}main{margin:0;padding:0;max-width:none}.sheet{padding:0;border:0;border-radius:0;margin:0;break-before:page}.sheet:first-child{break-before:auto}
h2{font-size:25pt;margin-bottom:8pt}.topic{font-size:9pt}h3{font-size:11.3pt;margin:8pt 0 3pt}p{margin:5pt 0}.sub,.small{font-size:8.6pt}.scroll{overflow:visible}th,td{font-size:8.9pt;padding:4pt 5pt}.note{padding:6pt 9pt}figure{margin:5pt 0}figure img{max-height:89mm;object-fit:contain}.chem img{max-height:95mm}.paired img{max-height:94mm}a{color:inherit}tr,figure{break-inside:avoid}}
"""


def build(root):
    root=Path(root)
    j=read(root/"results/comparison.json")
    pairs=read(root/"data/pair/summary.json")["replicas"]
    groups={"HH":"H/H","MeOH":"α-Me/ω-OH"}
    metrics={k:v["metrics"] for k,v in j["cases"].items()}
    pilot=j["archived_pilot"]["metrics"]
    def pm(entry, digits):
        return f'{tr(entry["mean"],digits)} ± {tr(entry["between_run_sd"],digits)}'
    def pmrow(title,key,digits=2):
        return [title,*(pm(metrics[k][key],digits) for k in groups)]
    duration=j["matched_duration_ps"]
    rawlink="https://github.com/burhanbeycan/burhanbeycan.github.io/tree/polymer-research-review-v0.2.0/polymer-research/polymer-solvent-md/hydration"
    nav='<nav><div>'+''.join(f'<a href="#{k}">{v}</a>' for k,v in [("model","Zincirler"),("size","Örnekleme"),("water","Hidratasyon"),("exchange","Su değişimi"),("pair","Su girişi"),("evidence","Yöntem ve veri")])+'</div></nav>'
    intro=f'''<p>Bir PEtOx trimerinin sudaki boyutunu ve karbonil gruplarının hidratasyonunu uç kimyasıyla ilişkilendirmek için iki sonlu zinciri aynı sıcaklık, basınç ve su sayısında örnekledik. Ayrı bir hesapta temas eden iki trimerin arasına su girişini atom kimliklerini koruyarak izledik.</p>
    {picture(root/"results/end_groups.png","Aynı PEtOx tekrarları ve iki farklı uçlandırmanın kimyasal yapıları.","chem")}
    {table(["Model","Temel soru","Koşullar"],[
    ["Tek zincir · H/H","Önceki 24 ps’lik gözlem, örnekleme uzatılınca nasıl değişiyor?",f"DP 3 · 500 su · üç koşu × {tr(duration,0)} ps"],
    ["Tek zincir · α-Me/ω-OH","Aynı DP’de iki uç birlikte değişince hangi bölgelerin su çevresi değişiyor?",f"DP 3 · 500 su · üç koşu × {tr(duration,0)} ps"],
    ["Temas çifti · H/H","İki zincirin yakın temasına su ne zaman ve hangi atomların arasından giriyor?",f"2 × DP 3 · 700 su · iki koşu × {tr(pairs[0]['duration_ps'],0)} ps"]])}
    <p class=note>α-Me/ω-OH seçimi literatürde bildirilen metil tosilat başlatması ve NaOH sonlandırmasına dayanır. Sıcaklık 298,15 K, hedef basınç 1 bar. İki uç modeli aynı anda değiştiğinden farklar uç çiftine aittir.</p>'''
    parts=[section("model","Polymer Research Studio · Açık su MD","Uç grupları ve suya açılan zincir teması",intro)]
    rgrows=[]
    for name,m in [("H/H · ilk 24 ps",pilot),("H/H · uzatılmış",metrics["HH"]),("α-Me/ω-OH · uzatılmış",metrics["MeOH"])]:
        r=m["rg_nm"]; ac=[v["ess"] for v in r["per_run_acf"] if v["ess"] is not None]
        rgrows.append([name," / ".join(tr(x,3) for x in r["per_run_mean"]),pm(r,3),tr(r["split_rank_rhat"],2),
                       (tr(min(ac),1)+"–"+tr(max(ac),1)) if ac else "—"])
    measured=[pmrow("Bütün zincir · Rg / nm","rg_nm",3),pmrow("Ortak 21 ağır atom · Rg / nm","core_rg_nm",3)]
    back=read(root/"results/backbone.json")["metrics"]
    for title,key in [("Ortak 9 omurga atomu · Rg / nm","backbone_rg"),("α-N ile son omurga C arası / nm","alpha_n_to_omega_backbone_c")]:
        measured.append([title,*[tr(back[c][key]["mean_nm"],3)+" ± "+tr(back[c][key]["between_lineage_sd_nm"],3) for c in ("H/H","Me/OH")]])
    body=f'''<p>Üç H/H uzatması, arşivlenmiş üç pilot koşunun son koordinat ve kutularından başladı. Kaydedilmemiş hızlar yeniden çekildi. 12 ps yeniden koşullandırmadan sonraki kısım analiz edildi. Karşılaştırma modeline aynı hazırlık süreleri uygulandı.</p>
        {picture(root/"results/convergence.png","Aynı fiziksel zamanda bağımsız koşularda zincir boyutu ve blok davranışı.")}
        {table(["Model","Koşu ortalamaları · nm","Rg · nm","R̂","ESS / koşu"],rgrows)}
        <p>Rg, bağlı zincir periyodik sınırlar boyunca bütünleştirildikten sonra kütle ağırlıklı hesaplandı. ± üç koşu ortalamasının örnek standart sapmasıdır. R̂, sıra-normalize bölünmüş ve katlanmış tanıların büyük olanıdır. ESS ardışık örneklerin özilintisini hesaba katan yaklaşık etkin örnek sayısıdır.</p>
        {table(["Boyut tanımı","H/H","α-Me/ω-OH"],measured)}
        <p class=note>H/H için R̂ ilk 24 ps’de {tr(pilot["rg_nm"]["split_rank_rhat"])}, uzatmada {tr(metrics["HH"]["rg_nm"]["split_rank_rhat"])} çıktı. Me/OH için {tr(metrics["MeOH"]["rg_nm"]["split_rank_rhat"])} bulundu. Üretim beş katına çıktığında da denge zincir boyutu doğrulanmadı. Hidratasyon ve zincir boyutu farklı hızlarda örnekleniyor.</p>'''
    parts.append(section("size","Tek zincir · Zaman ve örnekleme","Zincir boyutu koşular arasında tekrarlanıyor mu?",body))
    rows=[]
    for prefix,lab in [("alpha","α tarafı"),("middle","Orta birim"),("omega","ω tarafı")]:
        rows.append([lab,*(pm(metrics[c][prefix+"_neighbors"],2) for c in groups),
                     *(pm(metrics[c][prefix+"_water_hb"],2) for c in groups)])
    water=f'''<p>Karbonil oksijenleri α ucundan ω ucuna sıralandı. Her bölgede 0,35 nm içinde bulunan su oksijenlerini saydık. Hidrojen bağı için ayrıca su O–H···karbonil O açısının en az 150° olmasını aradık.</p>
        {picture(root/"results/hydration.png","Karbonillerin bölgeye ve uç kimyasına göre su komşuları ve yönlü hidrojen bağları.")}
        {table(["Karbonil","Komşu su · H/H","Komşu su · Me/OH","H bağı · H/H","H bağı · Me/OH"],rows)}
        {table(["Ortalama / bölge","H/H","α-Me/ω-OH"],[pmrow("Komşu su / karbonil","carbonyl_neighbors"),pmrow("Su → karbonil H bağı","water_to_carbonyl_hb"),pmrow("Uç donör → su H bağı¹","terminal_to_water_hb")])}
        <p class=sub>¹ H/H’de α-amid N–H, Me/OH’de ω-alkol O–H donördür. Kapanma ile üçüncül hale gelen α-amid N donörlüğünü kaybeder. Eklenen alkol donör veya akseptör olabilir. Kaydedilen su → ω-OH H bağı, ω-OH başına {pm(metrics["MeOH"]["water_to_terminal_OH_hb"],2)}. Tek su farklı karbonillere komşu olabilir. ω karbonilinde su vericili H bağı ortalaması {tr(metrics["HH"]["omega_water_hb"]["mean"])} iken uçlandırılmış modelde {tr(metrics["MeOH"]["omega_water_hb"]["mean"])} oldu. Alkol oksijeni ayrıca yeni bir hidratasyon bölgesi sunar.</p>'''
    parts.append(section("water","Tek zincir · Bölgesel etkileşim","Hangi karbonilin çevresi değişiyor?",water))
    retention=[]
    for site,label in [("alpha","α"),("middle","Orta"),("omega","ω")]:
        row=[label]
        for g in groups:
            item=next(x for x in j["cases"][g]["water_exchange"] if x["site"]==site and x["lag_ps"]==10)
            row.extend([tr(100*item["endpoint_mean"],1),tr(100*item["all_sampled_mean"],1)])
        retention.append(row)
    exchange=f'''<p>Karbonil çevresinde aynı anda iki su bulunması, aynı iki suyun orada kaldığı anlamına gelmez. Her zaman başlangıcında komşu olan suyun kimliğini kaydettik. 1, 2, 5 ve 10 ps sonra aynı karbonilin 0,35 nm kabuğunda bulunup bulunmadığını denetledik.</p>
        {picture(root/"results/water_exchange.png","Aynı suyun karbonil çevresinde yeniden görülmesi ve tüm örneklenmiş karelerde orada kalması.")}
        {table(["Karbonil","H/H · 10 ps sonra %","H/H · tüm karelerde %","Me/OH · 10 ps sonra %","Me/OH · tüm karelerde %"],retention)}
        <p class=note>“10 ps sonra” koşulu arada çıkıp yeniden giren suyu da sayar. “Tüm karelerde” koşulu aradaki her kayıtlı kareyi denetler. Kareler arasında 1 ps vardır. Bu yüzden sürekli mikroskobik ikamet veya H-bağı ömrü tayin edilmedi.</p>'''
    parts.append(section("exchange","Tek zincir · Su molekülünün kimliği","Aynı su karbonil çevresinde ne kadar kalıyor?",exchange))
    pairrows=[]
    for r in pairs:
        fir=r["windows"]["first_quarter"];end=r["windows"]["last_quarter"]
        pairrows.append([str(r["seed"]),tr(fir["contact_pairs"],1)+" → "+tr(end["contact_pairs"],1),
                         tr(fir["shared_waters"],1)+" → "+tr(end["shared_waters"],1),str(r["samples_with_zero_contact"]),
                         tr(r["interposition_first_ps"],0) if r["interposition_first_ps"] is not None else "Gözlenmedi"])
    both=f'''<p>İki H/H trimerini temas halinde yerleştirip suyu paketledik. Başlangıçtaki 8 ps koşullandırmadan sonra konum kısıtları kaldırıldı. 0 ps serbest hareketin başlangıcıdır. İki farklı başlangıç yerleşimi kullanıldı.</p>
        {picture(root/"results/water_entry.png","Başlangıçta temas eden iki ağır atomun arasına giren belirli bir suyun gerçek koordinatlardaki ölçümleri.","paired")}
        {table(["Tohum","Temas çifti¹","Ortak su¹","Temassız kare","İlk yerel su girişi · ps"],pairrows)}
        <p class=sub>¹ İlk ve son zaman çeyreğinin ortalamaları. Temas, iki zincirden birer ağır atomun ≤0,45 nm uzaklıkta olmasıdır. Ortak su, oksijeni her iki zincirin de birer ağır atomuna ≤0,35 nm uzaktaki sudur.</p>
        <p class=note>Yerel su girişi için başlangıçta temas eden belirli atom çifti açılmalı, aynı su O her iki atoma ≤0,35 nm yaklaşmalı ve zincir atomu–O–zincir atomu açısı ≥120° olmalıdır. Bu geometrinin ilk kaydedildiği zaman 1 ps örnekleme aralığıyla sınırlıdır. İki trimerlik sistem bir polimer katısının çözünme hızını ölçmez. Her iki koşuda son karede temas sürdü. İkinci koşu 44 ve 45 ps’de 0,45 nm ölçütüne göre temassızdı, 46 ps’de yeniden temas görüldü. 0,50 nm eşiğinde bu kısa kesinti de kayboldu.</p>
        <p><a class=button href="contact-viewer.html">İki koşunun 3B filmini aç</a></p>'''
    parts.append(section("pair","İki zincir · Temasın ıslanması","Su, açılan bir zincir temasına giriyor mu?",both))
    with (root/"results/amide_torsions.csv").open() as fp:
        torsions=list(csv.DictReader(fp))
    angle_rows=[]
    for name in ("HH","MeOH"):
        for site,lab in [("alpha","α"),("middle","Orta"),("omega","ω")]:
            group=sorted([r for r in torsions if r["chemistry"]==name and r["site"]==site],key=lambda r:int(r["seed"]))
            orientations=["0°" if float(r["frac_abs_phi_le_60"])==1 else "180°" if float(r["frac_abs_phi_ge_120"])==1 else "Karışık" for r in group]
            angle_rows.append([groups[name],lab,*orientations,sum(int(r["sampled_basin_changes"]) for r in group)])
    evidence=f'''<h3>Amid yönelimi ve başlangıç konformeri</h3>
        <p>O(karbonil)–C(karbonil)–N(amid)–C(ω yönündeki ilk metilen) burulmasını her bölgede ölçtük. Her satırda üç koşunun ziyaret ettiği bölgeler gösteriliyor. 0° sınıfı |φ| ≤60°, 180° sınıfı |φ| ≥120° olarak tanımlandı.</p>
        {table(["Zincir","Bölge","13031","27011","49033","Gözlenen geçiş"],angle_rows)}
        <p class=note>Orta amid H/H’de 0°, Me/OH’de 180° çevresinde kaldı. Her iki kimyada ω-amid yönelimi koşuya bağlıydı. Bu örneklemede zincir boyutu ve hidratasyon farkları, seçili amid yönelimlerinin etkisini de taşır. Sonraki sınama, her iki kimyayı da her iki orta-amid yöneliminden başlatmaktır.</p><h3>Moleküler model</h3><p>GAFF2 ve tüm sonlu moleküle hesaplanan AM1-BCC yükleri, rijit TIP3P su, OpenMM, periyodik NPT, PME, 0,8 nm gerçek uzay kesimi, 2 fs integrasyon adımı ve hidrojen içeren bağlara kısıt kullanıldı. GAFF2 aktarımları ve başlangıç koordinatları iki kimya için arşivlidir. Tüm koşular klasik MD ile yürütüldü.</p>
        <h3>Karşılaştırmanın kapsamı</h3><p>Üç bağımsız başlangıç yalnızca aynı kimyadaki koşular arasında tekrar oluşturur. H/H ve Me/OH fiziksel olarak eşleştirilmiş deney çiftleri değildir. Aynı su sayısı ve DP’de H/H zincirin kütle payı %3,217, α-Me/ω-OH zincirin kütle payı %3,528’dir. Uçların ayrı etkisini sınamak için α-Me/ω-H ve α-H/ω-OH çapraz kontrolleri gerekir. Daha uzun zincir ve daha büyük kutu ile duyarlılık sınaması, oligomerden polimere aktarımın sonraki adımlarıdır.</p>
        <h3>Veriye ulaşma</h3><p>İndirilebilir pakette tam atom koordinatları, periyodik kutular, yeniden başlatma bilgileri, kuvvet alanı paketleri, koşu ve dosya özetleri bulunur. Karşılaştırma CSV’leri her koşunun zaman serilerini, karbonil bazında ortalamaları, zaman bloklarını ve özilinti tanılarını içerir. Etkileşimli film aynı kaydedilmiş karelerde molekülleri ve sayısal izleri eşzamanlı değiştirir.</p>
        <ul><li><a href="{rawlink}/results">Sayısal sonuçlar ve grafikler</a></li><li><a href="{rawlink}/data">Kaynak topolojiler, parametreler ve yörüngeler</a></li><li><a href="{rawlink}/experiments">Koşu ve analiz betikleri</a></li></ul>
        <h3>Yöntem kaynakları</h3><p>PEtOx’un CROP uç kimyası · <a href="https://doi.org/10.1039/D2PY00437B">Polymer Chemistry, 2022</a>. Moleküler simülasyonlarda belirsizlik ve örnekleme kalitesi · <a href="https://doi.org/10.33011/livecoms.1.1.5067">Grossfield ve ark., 2018</a>. Sıra-normalize R̂ tanısı · <a href="https://doi.org/10.1214/20-BA1221">Vehtari ve ark., 2021</a>. Polimerlerin çözünme mekanizmaları · <a href="https://doi.org/10.1016/S0079-6700(03)00045-5">Miller-Chou ve Koenig, 2003</a>.</p>'''
    parts.append(section("evidence","Yöntem · Sınırlamalar · Veri","Hesap kaydı ve yeniden üretim",evidence))
    doc='<!doctype html><html lang="tr"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PEtOx · Uç grupları ve suya açılan temas</title><style>'+STYLE+'</style><header><div><span>Polymer Research Studio</span></div></header>'+nav+'<main>'+''.join(parts)+'</main></html>'
    (root/"comparison.html").write_text(doc,encoding="utf-8")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("root",type=Path)
    build(parser.parse_args().root)
