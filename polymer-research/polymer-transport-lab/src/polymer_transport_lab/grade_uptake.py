"""Blocked comparison of independently measured water uptake and PVA descriptors."""
import csv
import io
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.resources import files

import numpy as np

from .model import finite_number

MODE = "pva_grade_uptake"
DATA = files("polymer_transport_lab").joinpath("data")
SPEC = ("specimen,grade,md_origin,dh_fraction,Mn_g_mol,Xc_fraction,mixed_dyad,water_per_oh,"
        "water_per_oac,md_wref_g_g,md_T_K,md_protocol,nmr_source,md_source").split(",")
OBS = "specimen,study,aw,T_K,water_g_g,source".split(",")
LEVELS = ("Composition + conditions", "+ quantitative NMR dyads", "+ MD site hydration")
POLICY = ("grade", "study", "study_new_grade")


def csv_table(body, header, limit):
    if not isinstance(body, str) or len(body) > 500000:
        raise ValueError("Supply a CSV file as text, up to 500 kB.")
    reader = csv.DictReader(io.StringIO(body.strip()))
    if reader.fieldnames != header:
        raise ValueError("Expected CSV header: " + ",".join(header))
    rows = []
    for row in reader:
        if None in row or any(v is None for v in row.values()):
            raise ValueError("Each CSV row must contain the declared number of columns.")
        rows.append({k: v.strip() for k, v in row.items()})
        if len(rows) > limit:
            raise ValueError(f"At most {limit} records per table.")
    if not rows:
        raise ValueError("The CSV has no data rows.")
    return rows


def optional(row, key, lo, hi):
    return finite_number(row[key], key, lo, hi) if row[key] else None


@dataclass
class Specimen:
    name: str
    grade: str
    dh: float
    mass: float
    xc: float
    mixed: float | None
    hydration: tuple | None
    protocol: tuple | None

    def features(self, activity, temperature):
        # Per-repeat bulk faces only, equilibrated adsorption branch, no salt/additives.
        # Include both activity and the log activity odds for curvature in dilute/high-RH regimes.
        base = [activity, math.log(activity / (1-activity)), 1000 / temperature, self.dh,
                math.log(self.mass), self.xc]
        pair = [self.mixed, activity*self.mixed]
        oh, oac = self.hydration
        return base, base + pair, base + pair + [oh, oac, activity*oh, activity*oac]


def read_specimens(body):
    specimens, protocols, origins = {}, set(), set()
    for row in csv_table(body, SPEC, 300):
        name, grade = row["specimen"], row["grade"]
        if not name or not grade or name in specimens:
            raise ValueError("Specimen IDs must be nonempty and unique, with an explicit grade.")
        dh = finite_number(row["dh_fraction"], "OH repeat fraction", 0.000001, .999999)
        mass = finite_number(row["Mn_g_mol"], "Number-average molar mass", 1, 1e8)
        xc = finite_number(row["Xc_fraction"], "Dry PVA crystallinity fraction", 0, 1)
        mixed = optional(row, "mixed_dyad", 0, 1)
        if mixed is not None:
            if not row["nmr_source"] or mixed > 2*min(dh, 1-dh) + 1e-6:
                raise ValueError("Mixed-dyad fraction must fit composition and carry its NMR source.")
        hyd = [optional(row, key, 0, 100) for key in ("water_per_oh","water_per_oac")]
        protocol = None
        if any(v is not None for v in hyd):
            if None in hyd or not row["md_protocol"] or not row["md_source"]:
                raise ValueError("MD data need both site coordinations, a protocol and a source.")
            ref = finite_number(row["md_wref_g_g"], "Prespecified MD water loading g/g", 1e-8, 10)
            temp = finite_number(row["md_T_K"], "MD reference temperature K", 200, 450)
            if row["md_origin"] not in ("independent_reference","synthetic"):
                raise ValueError("MD reference loading must be fixed independently of uptake outcomes.")
            protocol = (ref, temp, row["md_protocol"])
            protocols.add(protocol)
            origins.add(row["md_origin"])
        specimens[name] = Specimen(name, grade, dh, mass, xc, mixed,
                                   tuple(hyd) if None not in hyd else None, protocol)
    if len(protocols) > 1 or len(origins) > 1:
        raise ValueError("Use a common independent reference loading, temperature and MD protocol.")
    return specimens, origins


def read_observations(body, specimens):
    obs, cells = [], set()
    for row in csv_table(body, OBS, 1000):
        name, study = row["specimen"], row["study"]
        if name not in specimens or not study or not row["source"]:
            raise ValueError("Every film measurement needs a registered specimen, study and original source.")
        aw = finite_number(row["aw"], "Interior water activity", 1e-6, 1-1e-6)
        t = finite_number(row["T_K"], "Vapour experiment temperature K", 200, 450)
        w = finite_number(row["water_g_g"], "Measured positive dry-basis uptake", 1e-10, 10)
        key = (name, aw, t)
        if key in cells:
            raise ValueError("A specimen/activity/temperature is duplicated. Average technical replicates and cite them.")
        cells.add(key)
        obs.append((specimens[name], study, aw, t, w))
    return obs


def ridge_predict(x, y, weights, x_test):
    """Weighted ridge, fixed penalty one, unpenalized intercept, train-only scaling."""
    mu = np.average(x, axis=0, weights=weights)
    scale = np.sqrt(np.average((x - mu)**2, axis=0, weights=weights))
    scale[scale < 1e-10] = 1
    train = (x-mu) / scale
    mean_y = np.average(y, weights=weights)
    wt = train.T * weights
    coef = np.linalg.solve(wt@train + np.eye(x.shape[1]), wt@(y-mean_y))
    return mean_y + (x_test-mu)/scale@coef


def compare(spec_text, obs_text, policy, origin):
    spec, origins = read_specimens(spec_text)
    obs = read_observations(obs_text, spec)
    if origin not in ("synthetic","user_supplied") or policy not in POLICY:
        raise ValueError("Choose a supported holdout and data origin.")
    if origin == "user_supplied" and ("synthetic" in origins or any(r[1].startswith("fixture-") for r in obs)):
        raise ValueError("Synthetic example records must retain the synthetic origin.")
    usable = [r for r in obs if r[0].mixed is not None and r[0].hydration is not None]
    groups = sorted({r[0].grade if policy == "grade" else r[1] for r in usable})
    allspec = {r[0].name for r in obs}
    complete = {r[0].name for r in usable}
    coverage = [[n,spec[n].grade, sum(r[0].name==n for r in obs),
                 "available" if spec[n].mixed is not None else "missing",
                 "available" if spec[n].hydration is not None else "missing",
                 "included" if n in complete else "excluded"] for n in sorted(allspec)]
    folded, predictions, skipped = [], [], []
    for group in groups:
        getgroup = (lambda r: r[0].grade) if policy == "grade" else (lambda r: r[1])
        test = [r for r in usable if getgroup(r)==group]
        test_names = {r[0].name for r in test}
        test_grades = {r[0].grade for r in test}
        train = [r for r in usable if getgroup(r)!=group and r[0].name not in test_names
                 and (policy!="study_new_grade" or r[0].grade not in test_grades)]
        if len({r[0].grade for r in train}) < 3 or (policy != "grade" and len({r[1] for r in train}) < 2):
            skipped.append(group)
            continue
        per_film = Counter(r[0].name for r in train)
        weights = np.array([1/per_film[r[0].name] for r in train])
        train_y = np.log([r[4] for r in train])
        obs_y = np.array([r[4] for r in test])
        for i, label in enumerate(LEVELS):
            x = np.array([r[0].features(r[2],r[3])[i] for r in train])
            xt = np.array([r[0].features(r[2],r[3])[i] for r in test])
            log_pred = ridge_predict(x,train_y,weights,xt)
            if (abs(log_pred)>700).any():
                raise ValueError("Holdout prediction exceeds numerical range. Inspect training/holdout coverage.")
            pred = np.exp(log_pred)
            per_specimen = defaultdict(list)
            for row, measured, value, feature in zip(test,obs_y,pred,xt):
                per_specimen[row[0].name].append(abs(value-measured))
                outside = ((feature<x.min(axis=0)-1e-9)|(feature>x.max(axis=0)+1e-9)).any()
                predictions.append([group,label,row[0].name,row[2],row[3],float(measured),float(value),
                                    "outside" if outside else "inside"])
            mae = float(np.mean([np.mean(e) for e in per_specimen.values()]))
            folded.append([group,label,len(per_film),len(per_specimen),mae])
    return coverage, folded, predictions, skipped


EVIDENCE = [
    ["Nakajima & Gocho 1978","Homogeneous and heterogeneous acetylation","Vapour-equilibrium isotherms",
     "Route and crystallinity co-vary; quantitative NMR dyads and MD coordinations absent from extracted observations."],
    ["Highly blocky AcPVOH 2026","Composition-matched random and blocky series","¹³C-NMR dyads",
     "The sequence study supplies no matched equilibrium-vapour/MD labels in the extracted data."],
    ["Chiessi et al. 2007","Chemically crosslinked PVA hydrogel at high hydration","Atomistic MD and QENS",
     "Network at 323 K; same-grade partially acetylated film and vapour equilibrium labels unavailable."],
]
REFS = [
    {"id":"pva_sequence_vapour","title":"Nakajima and Gocho (1978), equilibrium water-vapour sorption in differently acetylated polymers","url":"https://doi.org/10.1246/nikkashi.1978.1431"},
    {"id":"pva_matched_nmr","title":"Synthesis and Sequence Analysis of Highly Blocky Acetylated PVA Copolymers (2026), matched series and NMR","url":"https://doi.org/10.1021/acsapm.6c00317"},
    {"id":"pva_qens","title":"Water and polymer dynamics in chemically cross-linked PVA hydrogels (2007), atomistic MD and QENS","url":"https://doi.org/10.1021/jp0671143"},
]


def controls(app, select):
    vis = {"key":"mode","value":MODE}
    app["fields"][0]["options"].append({"value":MODE,"label":"PVA · does sequence or hydration improve equilibrium uptake?"})
    app["fields"].extend([
        select("grade_data","Dataset","papers",[("papers","Primary evidence inventory"),("csv","Two matched CSV tables"),("coordinates","Extract site hydration from coordinates")],
               "For known primary evidence, inspect the limits on cross-study matching.","PVA sequence and hydration",show_if=vis),
    ])
    table_vis = {"key":"grade_data","value":"csv"}
    for key, label, file in [("grade_spec","Specimen, NMR and MD reference records","pva_grade_spec.csv"),("grade_obs","Equilibrium film adsorption observations","pva_grade_obs.csv")]:
        app["fields"].append({"key":key,"label":label,"type":"textarea","unit":"CSV","default":DATA.joinpath(file).read_text(),
                              "group":"Matched records","help":"Edit or paste CSV, preserving the header and units.","show_if":table_vis})
    app["fields"].extend([
        select("grade_holdout","External split","grade",[("grade","Leave one complete grade out"),
               ("study","Leave one complete study out"),("study_new_grade","Leave study out and exclude its grades")],
               "A study holdout also purges any shared specimen; choose the joint split to reserve unseen grades.","Validation",show_if=table_vis),
        select("grade_origin","Data origin","synthetic",[("synthetic","Synthetic worked example"),("user_supplied","User-supplied records")],
               "Example coordinations and observations are generated placeholders.","Validation",show_if=table_vis),
    ])
    coord_vis={"key":"grade_data","value":"coordinates"}
    app["fields"].extend([
        {"key":"grade_sites","label":"Selected oxygen trajectories (CSV)","type":"textarea","unit":"nm","default":DATA.joinpath("pva_site_fixture.csv").read_text(),"help":"Enter one selected oxygen per frame and an orthorhombic cell. Roles: oh, oac (carbonyl O), water (oxygen).","group":"Hydration trajectories","show_if":coord_vis},
        {"key":"grade_cutoff","label":"Water-oxygen cutoff","type":"number","unit":"nm","default":0.35,"min":0.01,"max":1,"help":"Use and document the same radial definition across specimens.","group":"Hydration trajectories","show_if":coord_vis},
    ])
    app["presets"].extend([
        {"name":"PVA · evidence for sequence and uptake","description":"Original directions in vapour sorption and NMR availability across publications.","inputs":{"mode":MODE,"grade_data":"papers"}},
        {"name":"PVA · oxygen coordination across a periodic boundary","description":"Manually placed oxygen sites exercise periodic distances and equal-replica averaging.","inputs":{"mode":MODE,"grade_data":"coordinates"}},
        {"name":"PVA · withhold complete grades (synthetic)","description":"Compare three models on the same unseen grade with two CSVs of artificial observations.","inputs":{"mode":MODE,"grade_data":"csv","grade_holdout":"grade","grade_origin":"synthetic"}},
    ])
    app["references"].extend(REFS)
    app["models"].append({"name":"Grouped comparison of PVA film adsorption",
                        "equation":"ridge(log w): conditions and composition → add mixed dyads → add per-site hydration",
                        "assumptions":["Positive interior humidity adsorption measurements, g/g dry film","Three models share identical specimens and holdouts","MD coordinates are from one prespecified, outcome-independent reference loading and protocol"],"references":[r["id"] for r in REFS]})


def run(vals, metric, plot):
    out = dict(metrics=[],tables=[],plots=[],warnings=[],notes=[],material_updates={"parameters":[]},handoffs=[])
    if vals["grade_data"] == "coordinates":
        return hydration_view(vals,metric)
    if vals["grade_data"] == "papers":
        out["metrics"] = [metric("source_count","Primary sources audited",3,""),metric("matched","Complete merged cohorts confirmed",0,"")]
        out["tables"] = [dict(title="What can be linked from the reviewed sources",columns=["Source","Material","Observation","Missing link"],rows=EVIDENCE)]
        out["notes"] = ["Nakajima and Gocho report similar uptake in the two synthesis-route families at low RH, and greater uptake for the homogeneous series at high RH (Figs. 5–6). Their samples also differ in crystalline organization.",
                        "A route label does not determine a quantitative NMR mixed-dyad fraction. Hydration in a high-water crosslinked hydrogel does not supply a partially acetylated film's coordination at its reference water loading.",
                        "The three-model comparison awaits matched specimen observations. Select the synthetic preset to inspect the calculation."]
        return out
    if vals["grade_data"] != "csv":
        raise ValueError("Unrecognized grade dataset.")
    coverage, folded, pred, skipped = compare(vals["grade_spec"],vals["grade_obs"],vals["grade_holdout"],vals["grade_origin"])
    out["tables"].append(dict(title="Specimen coverage",columns=["Specimen","Grade","Observations","NMR","MD","Common cohort"],rows=coverage))
    groups = sorted({r[0] for r in folded})
    out["metrics"].append(metric("eligible_folds","Eligible held-out groups",len(groups),""))
    if folded:
        for i, label in enumerate(LEVELS):
            mae = float(np.mean([r[4] for r in folded if r[1]==label]))
            out["metrics"].append(metric(f"grade_mae_{i}",label+" · MAE",mae,"g/g","Mean of equal-film fold errors, then equal held-out groups."))
        out["tables"].extend([dict(title="Model comparison per holdout",columns=["Held out","Model","Train specimens","Test specimens","MAE g/g"],rows=folded),
                              dict(title="Blind predictions",columns=["Held out","Model","Specimen","Activity","K","Observed g/g","Predicted g/g","Feature range"],rows=pred)])
        out["plots"].append(plot("Held-out error by group", "Group number (see fold table)", "MAE (g/g dry film)",
            [dict(name=label,x=list(range(1,len(groups)+1)),y=[r[4] for r in folded if r[1]==label]) for label in LEVELS]))
    if not folded:
        out["warnings"].append("No eligible comparison. Supply three independent training grades per fold and, for a study holdout, two training studies with eligible records.")
    if skipped:
        out["warnings"].append("Ineligible held-out groups: "+", ".join(skipped))
    if vals["grade_origin"]=="synthetic":
        out["warnings"].append("Synthetic benchmark. All feature and target values in the bundled example are artificial; error changes exercise the comparison.")
    out["notes"] = ["Target is equilibrium water mass per dry film mass on the adsorption branch. Each model uses the same complete cases; all standardization and coefficients use training data only.",
                    "MD features count first-shell waters per site at a single reference water loading common to all grades. Using the measured equilibrium loading of a held-out specimen to prepare its simulation would leak the target.",
                    "Keep film processing, crystallinity assay, NMR assignments and MD force field/water model comparable. Held-out matter outside the training feature ranges is flagged.",
                    "No tracer or water self-diffusivity enters the fit, nor is a Fickian film diffusivity estimated from equilibrium uptake."]
    return out


def hydration_frames(body, cutoff):
    """Mean number of water oxygen sites within a sphere around each polymer O.

    Orthorhombic periodic cells. Atom rows may be wrapped or unwrapped.
    Scan each replicate independently so a longer run does not dominate its peers.
    """
    from scipy.spatial import cKDTree
    cutoff = finite_number(cutoff, "Coordination cutoff nm", 0.01, 1)
    coords = csv_table(body, ["replica","frame","atom","role","x_nm","y_nm","z_nm","Lx_nm","Ly_nm","Lz_nm"], 12000)
    frames = defaultdict(dict)
    boxes = {}
    for row in coords:
        rep, atom, role = (row[k] for k in ("replica","atom","role"))
        if not rep or not atom or role not in ("oh","oac","water"):
            raise ValueError("Each selected oxygen needs a replica, atom ID and role (oh/oac/water).")
        try:
            frame = int(row["frame"])
        except ValueError:
            raise ValueError("Frame indices must be nonnegative integers.") from None
        if frame < 0:
            raise ValueError("Frame indices must be nonnegative integers.")
        key = (rep, frame)
        box = tuple(finite_number(row[k],k,0.01,1e4) for k in ("Lx_nm","Ly_nm","Lz_nm"))
        if min(box) <= 2*cutoff or boxes.setdefault(key,box)!=box:
            raise ValueError("An orthorhombic box must be constant within the frame and exceed twice the cutoff on every axis.")
        if atom in frames[key]:
            raise ValueError("An atom ID appears twice in the same frame.")
        position = tuple(finite_number(row[k],k,-1e5,1e5) for k in ("x_nm","y_nm","z_nm"))
        frames[key][atom]=(role,position)
    topologies = {}
    values = defaultdict(list)
    for (rep, frame), atoms in sorted(frames.items()):
        topology = {key: value[0] for key,value in atoms.items()}
        if topologies.setdefault(rep,topology)!=topology:
            raise ValueError("Selected atom identities and roles must be stable within each replica.")
        roles = defaultdict(list)
        for role, pos in atoms.values():
            roles[role].append(pos)
        if set(roles)!={"oh","oac","water"}:
            raise ValueError("Each frame needs OH oxygen, acetate carbonyl oxygen and water oxygen sites.")
        box = np.array(boxes[(rep,frame)])
        tree = cKDTree(np.mod(roles["water"],box),boxsize=box)
        frame_avg = [np.mean(tree.query_ball_point(np.mod(roles[role],box),cutoff,return_length=True))
                     for role in ("oh","oac")]
        values[rep].append(frame_avg)
    per_rep = [[rep,len(rows),*map(float,np.mean(rows,axis=0))] for rep,rows in sorted(values.items())]
    return per_rep, list(map(float,np.mean([row[2:] for row in per_rep],axis=0)))


def hydration_view(vals, metric):
    per_rep, (oh,oac) = hydration_frames(vals["grade_sites"],vals["grade_cutoff"])
    warning = ("The bundled coordinates are an analytic periodic-cell fixture with manually placed atoms."
               if vals["grade_sites"].strip()==DATA.joinpath("pva_site_fixture.csv").read_text().strip() else
               "Interpretation depends on the supplied atom selections, force field, water model, sequence, equilibration and reference water loading.")
    return dict(metrics=[metric("water_per_oh","Water oxygens per hydroxyl O",oh,"per site"),
                         metric("water_per_oac","Water oxygens per acetate carbonyl O",oac,"per site")],
                tables=[dict(title="Equal-replica coordination estimates",columns=["Replica","Frames","Water / OH","Water / OAc carbonyl"],rows=per_rep)],
                plots=[],warnings=[warning],notes=[
                    "Distances use an orthorhombic periodic box and the shared cutoff in nm. Select a physically justified cutoff (for example the declared first RDF minimum) and use the same definition across grades.",
                    "Each number counts nearby water oxygens per polymer site. Repeat-frame means are averaged within each run and then equally over runs. These distances alone do not define hydrogen bonds.",
                    "Transfer these two means to the specimen CSV with a trajectory identifier, force field and independent reference loading/temperature. The reference MD state must be set before inspecting held-out film uptake."],
                material_updates={"parameters":[]},handoffs=[])
