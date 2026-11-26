"""Sample-held-out PVA sequence benchmarks and explicit aqueous PEtOx/EI charge states."""
from statistics import mean

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from .model import Oligomer, TEMPLATES, build_oligomer, molecular_view
from .sequences import fraction, pva_ensemble

PVA_URL = "https://doi.org/10.1021/acsapm.6c00317"
PVA_CITATION = "McCrary et al., ACS Applied Polymer Materials (2026), Table 1, DOI 10.1021/acsapm.6c00317."
PVAC_REPEAT = "[*:1]CC(OC(=O)C)[*:2]"
# Sample name, preparation route, nominal DA from sample ID (%); Table 1 η and run lengths.
# Missing η (bAcPVOH5, bAcPVOH25, PVOH) and the commercial comparison grade have no fit role.
SAMPLES = (
    ("rAcPVOH4",  "solution", 4, 0.72, 23.2, 1.5),
    ("rAcPVOH16", "solution", 16, 0.95, 6.1, 1.3),
    ("rAcPVOH23", "solution", 23, 0.95, 4.9, 1.3),
    ("bAcPVOH9",  "film", 9, 0.22, 51.0, 5.1),
    ("bAcPVOH12", "film", 12, 0.18, 45.9, 6.3),
    ("bAcPVOH15", "film", 15, 0.30, 22.3, 3.9),
)
COUNTERIONS = {"chloride": "[Cl-]", "bromide": "[Br-]", "acetate": "CC(=O)[O-]"}
UNITS = {"E": "N(C(=O)CC)CC", "I": "NCC", "P": "[NH2+]CC"}

def field(key, label, kind, default, mode, **kwargs):
    return dict(key=key, label=label, type=kind, default=default, group="Microstructure",
                show_if={"key":"chemistry_mode","value":mode}, **kwargs)

FIELDS = [
    field("matched_acetate", "Model acetate content (both routes)", "number", .15, "pva_benchmark",
          min=0.001,max=.999,step="any",unit="mole fraction",help="Common composition for the two model populations; experimental grades retain their own measured identities."),
    field("ei_pattern", "Chain from α to ω (E / I / P)", "text", "EEIPEEIP", "petox_hydrolysis",
          help="E = intact EtOx, I = neutral EI, P = protonated EI. Spaces are allowed. Up to 64 repeat units."),
    field("ei_anion", "Counterion for protonated EI", "select", "chloride","petox_hydrolysis",
          options=[{"value":k,"label":k.capitalize()} for k in COUNTERIONS]),
]
PRESETS = [
    {"name":"PVA · 2026 NMR and independent sample holdout",
     "description":"Evaluate sequence predictions, then compare two model populations at the same acetate content.",
     "inputs":{"chemistry_mode":"pva_benchmark","matched_acetate":.15,"ensemble_seed":23,"ensemble_dp":30,"ensemble_count":400}},
    {"name":"WBPU blend partner · PEtOx/EI protonation",
     "description":"Resolve hydrolysis position, protonation state, terminal connectivity and charge balance before a blend simulation.",
     "inputs":{"chemistry_mode":"petox_hydrolysis","ei_pattern":"EEIPEEIP","ei_anion":"chloride"}},
]
MODELS = [
    {"name":"Blind sequence-index comparison","equation":"η̂_i = mean(η_j for other grades j made by the same route); Bernoulli η = 1",
     "assumptions":["Six grades from a single paper; preparation route labelled independently of the NMR outcome",
                    "Separate controlled simulation uses literature η to constrain an ensemble at fixed composition"],"references":["pva_mccrary"]},
    {"name":"Explicit hydrolysis and protonation state","equation":"DH = (nI + nP)/DP ; qpolymer/e = nP ; nanion = nP",
     "assumptions":["Linear α-methyl / ω-OH chains","Protonation microstate supplied; the counterion balances charge"],
     "references":["petox_ends","wbpu_blends"]},
]
REFERENCES = [
    {"id":"pva_mccrary","title":"McCrary et al. (2026), NMR sequence analysis of randomly and blockily acetylated PVOH, Table 1", "url":PVA_URL},
    {"id":"wbpu_blends","title":"Beycan, WBPU / PEtOx and hydrolysed PEtOx fibre blends, MACRO 2024, abstract, p. 176",
     "url":"https://static1.squarespace.com/static/62c326b96a09b762cf8302ac/t/667e8cd4578728784d3601c9/1719569633844/AbstractBook-MACRO2024Final.pdf"}
]


def held_out_predictions(samples):
    rows = []
    for name, route, _, observed, _, _ in samples:
        training = [r[3] for r in samples if r[1] == route and r[0] != name]
        if not training:
            raise ValueError("Each sample requires another grade from the same preparation route.")
        predicted = mean(training)
        rows.append([name, route, len(training), observed, 1.0, predicted,
                     abs(1.0 - observed), abs(predicted - observed)])
    return rows


def pva_benchmark(values):
    x = fraction(values["matched_acetate"], "Common acetate fraction")
    if x in (0,1):
        raise ValueError("Use a copolymer composition strictly between zero and one.")
    folds = held_out_predictions(SAMPLES)
    obs_table = [[r[0],r[1],r[2],r[3],r[4],r[5]] for r in SAMPLES]
    models, views, ensembles = [], [], {}
    # The focus pair is nearest in nominal composition within the central region.
    for i in (1,5):
        name, _, _, eta, _, _ = SAMPLES[i]
        mix = 2 * (1-x) * x * eta
        target = [(1-x) - mix/2, mix, x - mix/2]
        ensemble = pva_ensemble("pva_nmr",None,target,values["ensemble_dp"],values["ensemble_count"],values["ensemble_seed"])
        seq = ensemble["examples"][0]
        chain = build_oligomer(TEMPLATES["pva"]["smiles"], seq, PVAC_REPEAT)
        view = molecular_view(chain)
        view.update(title=f"{name} η · simulated chain", comparison_group="pva_sequence_route",
                    caption=f"Same ensemble acetate fraction {x:.3f}; chain 1 is a finite seeded draw. A = OH, B = OAc. Sequence: {seq}. H/H model caps; tacticity unspecified.")
        views.append(view)
        models.append([name, x, eta, mix, ensemble["observed_dyads"][1],ensemble["sample_oh_fraction"],seq])
        ensembles[name] = ensemble
    route_error = mean(row[7] for row in folds)
    random_error = mean(row[6] for row in folds)
    return {"metrics":[{"key":"blind_route_mae","label":"Route model · held-out η MAE","unit":"η units","value":route_error},
                       {"key":"bernoulli_mae","label":"Bernoulli model · held-out η MAE","unit":"η units","value":random_error},
                       {"key":"nmr_grades","label":"Eligible experimental grades","unit":"grades","value":len(folds)}],
            "molecular_views":views,
            "plots":[{"title":"Equal-composition ensembles · OH count distribution","x_label":"OH repeats per chain",
                      "y_label":"Simulated chains","series":[{"name":n,"x":list(range(int(values["ensemble_dp"])+1)),
                      "y":[dict(e["histogram"]).get(k,0) for k in range(int(values["ensemble_dp"])+1)]} for n,e in ensembles.items()]}],
            "tables":[{"title":"Measured NMR sequence statistics · Table 1","columns":["Grade","Route","Nominal acetate (%)","Measured η","OH run (units)","Ac run (units)"],"rows":obs_table},
                      {"title":"One entire grade held out · η from route alone","columns":["Test grade","Preparation","Training count","Measured η","Bernoulli η","Route η","|error| Bernoulli","|error| route"],"rows":folds},
                      {"title":"Two model populations at controlled composition","columns":["η source grade","Fixed OAc fraction","η constraint","Target mixed fraction","Sampled mixed fraction","Sampled OH fraction","First sequence"],"rows":models}],
            "pva_benchmark":{"common_acetate":x,"sample_holdouts":folds,"models":models},
            "warnings":["Sample holdout within one publication. Repeated grades, laboratories or independent preparation routes are not external validation.",
                        "The common-composition ensembles are constrained by the measured η of each chosen grade. Their realized dyads are a sampling check. Each experimental grade retains its own nominal acetate percentage."],
            "notes":[PVA_CITATION,"Sample numbers denote nominal acetylation from the published sample IDs. Reported η and run lengths come directly from Table 1 and are not reconstructed from rounded nominal percentages.",
                     "Prediction consumes only the preparation route and training η. The held-out η, run lengths, dyads and sample label do not enter the fold fit.",
                     "Rows without measurable η are excluded; the commercial Poval raw material is outside the two synthesis routes.",
                     "Water coordination, crystallinity and film transport require independent observables; these NMR data constrain sequence. Sequence graphs display the first of the generated finite chains."],
            "handoffs":[],"material_updates":{}}


def _chain(pattern):
    mol = Chem.RWMol(Chem.MolFromSmiles("C"))
    mol.GetAtomWithIdx(0).SetProp("end_role","alpha")
    tail = 0
    for unit_index, symbol in enumerate(pattern):
        fragment = Chem.MolFromSmiles(UNITS[symbol])
        first = mol.GetNumAtoms()
        for atom in fragment.GetAtoms():
            atom.SetIntProp("unit_index",unit_index)
        mol.InsertMol(fragment)
        mol.AddBond(tail,first,Chem.BondType.SINGLE)
        if unit_index == 0:
            mol.GetAtomWithIdx(first).SetProp("end_role","alpha")
            mol.GetBondBetweenAtoms(tail,first).SetProp("end_role","alpha")
        tail = mol.GetNumAtoms()-1
    o = Chem.Atom("O")
    o.SetProp("end_role","omega")
    omega = mol.AddAtom(o)
    mol.AddBond(tail,omega,Chem.BondType.SINGLE)
    mol.GetAtomWithIdx(tail).SetProp("end_role","omega")
    mol.GetBondBetweenAtoms(tail,omega).SetProp("end_role","omega")
    result = mol.GetMol()
    Chem.SanitizeMol(result)
    if Chem.GetFormalCharge(result) != pattern.count("P"):
        raise ValueError("Proton count and covalent graph charge differ.")
    return result


def hydrolysis(values):
    raw = values["ei_pattern"]
    if not isinstance(raw,str) or len(raw)>1024:
        raise ValueError("Enter an explicit 1–64 unit E/I/P sequence.")
    pattern = "".join(raw.split()).upper()
    if not 1 <= len(pattern) <= 64 or set(pattern)-set(UNITS):
        raise ValueError("Enter 1–64 E (EtOx), I (neutral EI) or P (protonated EI) symbols.")
    anion = values["ei_anion"]
    if not isinstance(anion,str) or anion not in COUNTERIONS:
        raise ValueError("Choose chloride, bromide or acetate.")
    ion = Chem.MolFromSmiles(COUNTERIONS[anion])
    dp, n_ei, q = len(pattern), len(pattern)-pattern.count("E"), pattern.count("P")
    neutral = pattern.replace("P","I")
    states = [("Intact PEtOx", "E"*dp)]
    if neutral != states[-1][1]: states.append(("Neutral hydrolysed chain",neutral))
    if pattern != states[-1][1]: states.append(("Selected charged microstate",pattern))
    views, ledger = [], []
    for title, seq in states:
        mol = _chain(seq)
        view = molecular_view(Oligomer(mol, seq, {}))
        view.update(title=title,comparison_group="hydrolysis_microstates",
                    caption="α-Me at the left end of the input sequence; ω-OH at the right. Repeat numbers follow input order. 2D connectivity; counterions are counted separately.")
        views.append(view)
        ledger.append([title, seq, rdMolDescriptors.CalcMolFormula(mol), Descriptors.MolWt(mol), Chem.GetFormalCharge(mol), view["smiles"]])
    charge_balance = q + q*Chem.GetFormalCharge(ion)
    if q:
        ionview = molecular_view(Oligomer(ion,"",{}))
        ionview.update(title=f"Counterion · {anion} × {q}",caption="One ion is drawn. Multiplicity follows the polymer charge. No ion-contact distance is specified.")
        views.append(ionview)
    return {"metrics":[{"key":"hydrolysis","label":"Hydrolysed repeats","value":n_ei/dp,"unit":"mole fraction"},
                       {"key":"protonation","label":"Protonated fraction of EI sites","value":q/n_ei if n_ei else None,"unit":"fraction"},
                       {"key":"polymer_charge","label":"Covalent-chain charge","value":q,"unit":"e"},
                       {"key":"salt_charge","label":"Chain + counterions","value":charge_balance,"unit":"e"}],
            "plots":[], "molecular_views":views,
            "tables":[{"title":"Matched-DP scenarios","columns":["Model","Sequence α→ω","Covalent formula","Average M (g/mol)","Charge (e)","Canonical SMILES"],"rows":ledger},
                      {"title":"Chemical sites for a blend comparison","columns":["Site","Count","Role under the selected microstate"],"rows":[
                       ["EtOx amide carbonyl",dp-n_ei,"H-bond acceptor"],["Neutral EI –NH–",n_ei-q,"Donor and acceptor"],
                       ["Protonated EI –NH₂⁺–",q,"Cationic donor"],["Terminal OH",1,"Terminal alcohol donor and acceptor"],
                       [anion,q,"Monovalent counterion"],["Salt formula unit",1, ledger[-1][5]+("."+COUNTERIONS[anion])*q]]}],
            "hydrolysis":{"sequence":pattern,"dp":dp,"intact_etox":dp-n_ei,"neutral_ei":n_ei-q,"protonated_ei":q,
                          "polymer_charge_e":q,"counterion":anion,"anion_count":q,"formal_total_charge_e":charge_balance},
            "notes":["I positions remove the N-propionyl side group from the EtOx backbone. The α-methyl and ω-hydroxyl termini are fixed in all graphs. No backbone cleavage is assumed.",
                     "These are supplied microstates. Sequence, absolute hydrolysis degree, pH-dependent protonation and anion identity require sample-specific evidence.",
                     "For WBPU blending, compare equal DP, blend mass fraction, ionic strength and water content, with care to account for the polymer/salt mass basis. Record counterions and replicate configurations before MD.",
                     "Morgan et al. reported terminal-group effects for their oligomers; the explicit site pattern is a selected scenario."],
            "warnings":["Stoichiometric counterions balance formal charge. The drawing does not assign ion pairing, 3D geometry, interaction energy or compatibility."],
            "handoffs":[],"material_updates":{}}
