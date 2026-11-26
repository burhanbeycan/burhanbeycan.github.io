"""PEO/LiTFSI benchmark and trajectory-analysis scenarios for the Studio."""
import json

from .displacements import displacement_result
from .shell import metric, plot, result, summarize_shell, table


DOI = "10.1021/acsmacrolett.4c00368"
# Printed values, 300 K, first shell through 2.7 Å, paragraph for Figure 3.
COORDINATION = {
    10: [(3.03, 2.37), (4.20, 1.13)],
    18: [(4.53, 1.21), (5.32, 0.10)],
}
MODELS = ("Original OPLS", "Scaled salt charges (0.79) + refitted Li–O LJ")


def published_shell(eo_li):
    if str(eo_li) not in ("10","18"):
        raise ValueError("The published coordination comparison has EO:Li=10 or 18.")
    ratio = int(eo_li)
    (ether0, tfsi0), (ether1, tfsi1) = COORDINATION[ratio]
    rows = [[label, ether, tfsi, ether+tfsi, ether/(ether+tfsi)]
            for label,(ether,tfsi) in zip(MODELS,COORDINATION[ratio])]
    traces=[]
    for donor, index in (("PEO ether oxygen",0),("TFSI sulfonyl oxygen",1)):
        for model, label in enumerate(MODELS):
            traces.append(dict(name=f"{donor} · {label}",x=[10,18],y=[COORDINATION[r][model][index] for r in (10,18)]))
    return result([table(f"Reported MD shell at EO:Li={ratio}, T=300 K, cutoff 0.27 nm",
                         ["Force field","Ether O/Li","TFSI O/Li","Total O/Li","Ether fraction"],rows)],
                  [plot("Two reported compositions (line joins the tabulated points)", "EO:Li", "O atoms per Li", traces)],
                  [f"{DOI}, results accompanying Figure 3 and model parameters in Table 2. Coordination numbers were integrated from simulated radial distribution functions. The summed values use rounded reported components.",
                   "The modified model changed both ionic partial charges and Li–O Lennard-Jones parameters. The coordination shift is a joint model sensitivity. The underlying article reused a DFT-based parametrization; this benchmark needs no new quantum calculations.",
                   "Only aggregate coordination numbers were transcribed. A contact residence time and replicate uncertainty require trajectories and cannot be recovered from mean coordination."],
                  metrics=[metric("ether_delta","Δ ether coordination",ether1-ether0,"O/Li"),
                           metric("tfsi_delta","Δ TFSI coordination",tfsi1-tfsi0,"O/Li"),
                           metric("total_delta","Δ total shell coordination",(ether1+tfsi1)-(ether0+tfsi0),"O/Li")])


def contact_example():
    state="synthetic-shell"
    runs=[]
    for model, offsets in (("Fixture model A",(0,1)),("Fixture model B",(1,2))):
        for i,offset in enumerate(offsets):
            frames=[]
            for t in range(9):
                ether_refs = [0,1] if model.endswith("A") else [0,1,2]
                frames.append({"ether":[[0,o] for o in ether_refs if (t+offset+o)%5!=0],
                               "tfsi":[[0,0]] if (t+offset)%4<2 else []})
            runs.append(dict(state_id=state,model=model,replica=f"fixture-{i}",force_field_source=f"{model}, invented for algorithm testing",trajectory_source="Deterministic constructed shell events",frames=frames))
    return dict(state_id=state,temperature_K=300,n_li=1,n_ether_o=10,dt_ps=2,cutoff_nm=.27,max_lag=4,runs=runs)


def moment_example():
    s="synthetic-linear-moments"
    return dict(state_id=s,temperature_K=300,volume_nm3=20,n_li=4,frame="polymer_com",
                transport_charge="formal_1e", fit_ps=[10,80],
                runs=[dict(state_id=s,model="Constructed diffusive limit",replica="analytic-1",
                           force_field_source="No potential: analytic regression fixture",
                           trajectory_source="Constructed moments with exact linear slopes",
                           moments=[[t,.00006*t,.00012*t,.00036*t] for t in range(10,90,10)])])


CONTACT = json.dumps(contact_example(),separators=(",",":"))
MOMENTS = json.dumps(moment_example(),separators=(",",":"))
OPTIONS = [
    ("md_reference", "PEO/LiTFSI · published force-field comparison"),
    ("md_shell", "Li coordination and residence · contact time series"),
    ("md_conductivity", "3D ion displacement moments · conductivity"),
]
FIELDS = [
    {"key":"md_eo_li","type":"select","label":"Published EO:Li","default":"10",
     "options":[{"value":v,"label":f"{v}:1"} for v in ("10","18")],"group":"MD evidence","show_if":{"key":"model","value":"md_reference"}},
    {"key":"md_contact_record","type":"textarea","label":"Complete Li–oxygen contact time series (JSON)","default":CONTACT,
     "help":"Consecutive uniformly sampled frames, stable atom IDs, minimum-image cutoff in nm. Each model/replica has its own provenance. The initial data are synthetic.",
     "group":"MD inputs","show_if":{"key":"model","value":"md_shell"}},
    {"key":"md_moment_record","type":"textarea","label":"Polymer-frame lagged displacement moments (JSON)","default":MOMENTS,
     "help":"Enter lag in ps, per-Li MSD, per-TFSI MSD and extensive squared signed ionic displacement (all moments in nm²). Include at least six lag values in the selected fit window.",
     "group":"MD inputs","show_if":{"key":"model","value":"md_conductivity"}},
]
PRESETS = [
    {"name":"PEO shell · original versus scaled force field","description":"Published MD: ether and anion coordination of Li, 300 K. Select EO:Li=10 or 18.","inputs":{"model":"md_reference","md_eo_li":"10"}},
    {"name":"Continuous Li–oxygen contacts · synthetic fixture","description":"Four constructed contact series illustrate exchange and between-replica spread.","inputs":{"model":"md_shell","md_contact_record":CONTACT}},
    {"name":"Collective charge motion · analytic fixture","description":"Known slopes in an artificial neutral system give σ/σNE = 0.5. Replace with trajectory moments.","inputs":{"model":"md_conductivity","md_moment_record":MOMENTS}},
]


def md_result(p):
    if p["model"]=="md_reference":
        return published_shell(p["md_eo_li"])
    if p["model"]=="md_shell":
        out=summarize_shell(p["md_contact_record"])
        if p["md_contact_record"]==CONTACT:
            out["warnings"].append("These are constructed event sequences for a software check. They are not MD observations.")
    else:
        out=displacement_result(p["md_moment_record"])
        if p["md_moment_record"]==MOMENTS:
            out["warnings"].append("The input contains exact analytic moments, not a simulation or conductivity measurement.")
    return out
