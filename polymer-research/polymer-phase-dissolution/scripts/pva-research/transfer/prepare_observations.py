"""Build the measured aqueous apparent-viscosity stress test. No model-labelled data."""
import argparse
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
PVA = "10.3390/polym15071743"
HPAM = "10.1007/s13202-017-0375-6"
CMC = "10.35812/cellulosechemtechnol.2020.54.24"
ALG = "10.1111/j.1365-2621.2009.02088.x"
FIELDS = ["id", "material", "group", "grade", "curve", "T_K", "shear_s_inv",
          "conc_value", "conc_basis", "NaCl_g_L", "carboxy_sites_per_repeat",
          "eta_cP", "replicates", "source_doi", "source_locator"]
# Laribi et al., Table 4; measured CMC, excluding the CMH cellulose/hemicellulose mixtures.
CMC_STRESS_TEST = [
    (170, 46.19, 274.46, 1201), (200, 49.53, 269.46, 1108),
    (300, 46.89, 240.38, 903.03), (500, 43.16, 201.42, 694.57),
    (700, 45.67, 174.04, 585.02), (1000, 36.66, 138.05, 476.93)]
# Table 1 (first observed apparent-viscosity column, Pa s), not Cross fits.
ALG_STRESS_TEST = [(1, [0.116, 0.085, 0.057]), (2, [0.396, 0.303, 0.225]),
                   (3, [1.225, 0.950, 0.672])]

def prepare(pva_path, hpam_path):
    rows = []
    def append(mat, broad, grade, temp, shear, c, basis, salt, carboxy, eta, n, doi, src):
        curve = f"{mat}/{grade}/{c:g}{basis}/{salt}gL/{temp:.2f}K"
        rows.append(dict(zip(FIELDS, [f"{mat}-{len(rows)+1:04d}", mat, broad, grade,
                curve, f"{temp:.2f}", shear, c, basis, salt, carboxy, eta, n, doi, src])))
    with open(pva_path) as f:
        for r in csv.DictReader(f):
            c, tc, shear, stress = (float(r[k]) for k in
                ("wt_percent", "temperature_C", "shear_rate_s_inv", "shear_stress_Pa"))
            append("PVA","vinyl_alcohol","Merck-141351",273.15+tc,shear,c,"wt_percent",
                   0,0,1000*stress/shear,1,PVA,"Table 2; reported stress/rate")
    # Articled ppm basis preserved. Multiple trials at identical conditions contribute
    # a single target (arithmetic mean on physical, not logarithmic viscosity scale).
    grouped = {}
    with open(hpam_path) as f:
        for r in csv.DictReader(f):
            k=tuple(float(r[m]) for m in ("T_C","shear_s","c_ppm","salt_g_L"))
            grouped.setdefault(k,[]).append(float(r["eta_mPa_s"]))
    for (tc,shear,c,salt),v in grouped.items():
        append("HPAM","acrylamide","SNF-05mol-percent",tc+273.15,shear,c,"ppm_reported",
               salt,.05,sum(v)/len(v),len(v),HPAM,"Table 1+2; 10^log10(eta), replicate mean")
    for shear,*viscosities in CMC_STRESS_TEST:
        for ds,eta in zip([.43,.88,1.23],viscosities):
            append("CMC","polysaccharides",f"CMC-DS{ds:g}",298.15,shear,4,
                   "g_per_100mL",0,ds,eta,1,CMC,"Table 4 CMC; room temperature (25 ± 1 °C in Methods)")
    for c,conditions in ALG_STRESS_TEST:
        for tc,eta in zip([25,40,55],conditions):
            append("ALG","polysaccharides","sodium-alginate",tc+273.15,100,c,
                   "wt_percent",0,1,1000*eta,1,ALG,"Table 1 apparent viscosity at 100/s")
    return rows

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pva",type=Path,default=HERE/"sources/carrasco_table2.csv")
    parser.add_argument("--hpam",type=Path,default=HERE/"sources/akbari_tables1_2.csv")
    parser.add_argument("--out",type=Path,default=HERE/"transfer_observations.csv")
    a=parser.parse_args()
    with a.out.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS); writer.writeheader()
        writer.writerows(prepare(a.pva,a.hpam))
