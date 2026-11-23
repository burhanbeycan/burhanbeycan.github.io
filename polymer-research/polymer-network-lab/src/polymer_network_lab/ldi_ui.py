"""Polyurethane chemistry, trajectory and validation views for the Studio GUI."""
import json

from .ldi import BEYCAN, HAO, CLASSES, DEMO_SERIES, DEMO_TRAJECTORY, FIXTURE_NOTE, SERIES_NOTE
from .ldi import contact_analysis, contact_folds, feed_scenario
from .wbpu_ui import _choice, _metric, _table, _text


def _number(key, label, mode, default, unit, help_text):
    return {"key":key,"type":"number","default":default,"label":label,"unit":unit,"step":"any","min":0,"max":{"ldi_angle":180,"ldi_cutoff":0.6}.get(key,1e8),
            "show_if":{"key":"network_mode","value":mode},"group":"LDI polyurethane","help":help_text}


def extend_ldi(app):
    feed="ldi_feed"; contacts="ldi_contacts"; folds="ldi_folds"
    app["fields"][0]["options"].extend([
        {"value":feed,"label":"LDI–PEG/PCL feed and linkage bounds"},
        {"value":contacts,"label":"Urethane N–H contacts by oxygen chemistry"},
        {"value":folds,"label":"Compare recipe and MD descriptors against FTIR"}])
    app["fields"].extend([
        *[_number(key,name,feed,value,"mmol","Amount of molecules on the same dry-component basis.") for key,name,value in [
            ("ldi_peg","PEG diol",.05),("ldi_pcl","PCL diol",1),("ldi_isocyanate","LDI molecules",1.65),
            ("ldi_ext","Added extender",.5),("ldi_lys","L-lysine in aqueous phase",.2)]],
        _choice("ldi_extender","Extender",[("BDO","1,4-Butanediol"),("PDO","1,3-Propanediol"),("glycerol","Glycerol")],
                "BDO","LDI polyurethane","Glycerol supplies three OH groups per intact molecule.",feed),
        _text("ldi_residual","Residual NCO at dispersion (mmol), optional","",feed,"LDI polyurethane",
              "Use measured total mmol on this batch before aqueous L-lysine is added. An empty field keeps the experimental upper bound unspecified."),
        _text("ldi_recipe_origin","Molar-amount source","Illustrative input amounts; Beycan et al. establishes the chemistry only.",
              feed,"LDI polyurethane","Identify your batch, LDI ester identity, feed assays and reaction stages."),
        _text("ldi_traj","Typed trajectory (JSON)",DEMO_TRAJECTORY,contacts,"Trajectory",
              "Coordinates and orthorhombic box in nm; use one chain ID per molecule and validated urethane N,H pairs. See docs/LDI.md for GRO trajectories."),
        _number("ldi_cutoff","N···O distance cutoff",contacts,.35,"nm","Minimum-image heavy-atom distance."),
        _number("ldi_angle","Minimum N–H···O angle",contacts,150,"degrees","180° is linear; explicit hydrogens required."),
        _text("ldi_traj_source","Source and model state",FIXTURE_NOTE,contacts,"Trajectory","For simulated data record the trajectory, force field, chain assignment, hydration, temperature, equilibration and sampling interval."),
        _text("ldi_ftirs","Formulation-level replicates (CSV)",DEMO_SERIES,folds,"Matched data","Supply independent spectra and trajectory IDs, identical FTIR band model, drying, conditioning and morphology."),
        _text("ldi_ftir_source","Origin of spectra and simulations",SERIES_NOTE,folds,"Matched data","Identify study, primary assay, trajectory protocols and lab records.")])
    app["presets"].extend([
        {"name":"LDI / PEG / PCL · feed and branching","description":"Beycan formulation chemistry with illustrative molar inputs; vary the extender.","inputs":{"network_mode":feed}},
        {"name":"LDI / PEG / PCL · NH contact geometry","description":"Periodic geometric control with named carbonyl and ether oxygen sites.","inputs":{"network_mode":contacts}},
        {"name":"LDI / PEG / PCL · recipe vs MD / FTIR","description":"Synthetic eight-formulation check of grouped prediction; replace with matched film data.","inputs":{"network_mode":folds}}])
    app["references"].extend([{"id":"beycan2024ldi","title":"Beycan et al. — LDI/PEG/PCL waterborne polyurethane electrospinning","url":"https://doi.org/"+BEYCAN,"year":2024},
                              {"id":"hao2016ldi","title":"Hao et al. — Lysine-based waterborne polyurethane","url":"https://doi.org/"+HAO,"year":2016}])
    app["models"].extend([
        {"name":"LDI feed and aqueous lysine bounds","equation":"n(OH)=2(nPEG+nPCL)+fext next; urea from lysine ≤ min(residual NCO, 2 nlys)",
         "assumptions":["Nominal functionality; independently measured NCO is required for the bound."],"references":["beycan2024ldi","hao2016ldi"]},
        {"name":"Interchain urethane hydrogen bonds","equation":"E_k = (contacts_k / all contacts) / (eligible pairs_k / all eligible pairs)",
         "assumptions":["Site-inventory normalization at fixed topology; orthorhombic minimum image and explicit N–H donors."],"references":[]},
        {"name":"Recipe vs recipe plus contact benchmark","equation":"Ridge objective = Σ (y − ŷ)² + Σ β²",
         "assumptions":["Entire formulations held out. Training-fold feature scaling, fixed regularization strength, unpenalized intercept."],"references":[]}])
    app["tagline"]="LDI polyurethane chemistry, molecular contacts, film evidence and network swelling"


def calculate_ldi(v):
    result={"metrics":[],"tables":[],"plots":[],"warnings":[],"notes":[],"material_updates":{},"handoffs":[]}
    out=result["metrics"]
    if v["network_mode"]=="ldi_feed":
        if not v["ldi_recipe_origin"].strip():
            raise ValueError("Identify the formulation and source for entered amounts.")
        r=feed_scenario(v["ldi_peg"],v["ldi_pcl"],v["ldi_isocyanate"],v["ldi_ext"],v["ldi_extender"],
                        v["ldi_residual"],v["ldi_lys"])
        out.extend([_metric("ncooh","Nominal NCO/OH index",r["index"],"%","Before the aqueous stage."),
                    _metric("oh","Dry OH inventory",r["oh"],"mmol","Nominal site count."),
                    _metric("triol","Potential three-arm molecules",r["triol"],"mmol","Glycerol molecules, functionality three before reaction.")])
        if r["urea_upper"] is not None:
            out.append(_metric("upper","Maximum lysine–NCO additions",r["urea_upper"],"mmol","Upper bound from the supplied residual NCO and two amine sites per lysine."))
        result["tables"]=[_table("Staged feed",["Component","Molecule amount (mmol)","Reactive sites (mmol)","Site type"],r["rows"])]
        result["notes"]=[v["ldi_recipe_origin"],"Beycan et al.: PEG/PCL, LDI then aqueous L-lysine. Source abstract identifies BDO, PDO and glycerol variants; project input amounts are illustrative. Specify the actual lysine diisocyanate ester from its lot or source before atomistic modeling.",
                         "The diols and diisocyanate provide a linear-chain candidate. Three reacted OH groups on one glycerol permit branching; neither gelation nor the number of reacted sites follows from feed alone. Pre-dispersion chain ends may carry NCO or OH.",
                         "At the aqueous stage, lysine amines and water can compete for residual isocyanate. Urea via lysine requires NCO surviving to dispersion. Water consumption, protonation, accessibility and side reactions affect actual connectivity. Record conversion or residual NCO assay and gel/sol evidence."]
    elif v["network_mode"]=="ldi_contacts":
        source=v["ldi_traj_source"].strip()
        if not source or (v["ldi_traj"].strip()!=DEMO_TRAJECTORY and source==FIXTURE_NOTE):
            raise ValueError("Identify the origin and conditions of the changed trajectory.")
        meta=json.loads(v["ldi_traj"])
        if not isinstance(meta,dict) or not isinstance(meta.get("frames"),list):
            raise ValueError("Supply atom assignments and frames in the documented JSON format.")
        if not isinstance(meta.get("atom_chain"),list):
            raise ValueError("Supply one chain ID per atom.")
        if len(meta["frames"])>1000 or len(meta["atom_chain"])>1000:
            raise ValueError("Use the streaming GRO analyzer for over 1000 atoms or frames.")
        r=contact_analysis(meta,meta["frames"],v["ldi_cutoff"],v["ldi_angle"])
        out.extend([_metric("frames","Frames",r["frames"],"","Equally spaced input."),
                    _metric("contacts","Qualifying triplet observations",r["total"],"","May include bifurcated donors."),
                    _metric("nh","Urethane N–H donors",r["donors"],"","User-verified bonded pairs.")])
        if r["lag_survival"] is not None:
            out.append(_metric("lag","One-frame pair retention",r["lag_survival"],"fraction",f"Contact present again {r['stride']:g} ps later."))
        result["tables"]=[_table("Interchain N–H···O contacts",["Oxygen group","Typed oxygen sites","Eligible N···O pairs","Triplet observations","Contacts / donor / frame","Inventory share","Enrichment"],r["rows"])]
        result["tables"].append(_table("Saved-frame pair persistence",["Oxygen group","One-frame retention","Longest observed span (ps)","Episodes bracketed by absence","Mean bracketed span (ps)"],r["persistence"]))
        x=[row[0] for row in r["by_frame"]]
        result["plots"]=[{"title":"Qualifying interchain contacts by oxygen group","x_label":"Time (ps)","y_label":"Triplets / frame","series":[
            {"name":k,"x":x,"y":[row[1+CLASSES.index(k)] for row in r["by_frame"]]} for k in [g[0] for g in r["rows"]]]}]
        result["notes"]=[source,"Enrichment 1 corresponds to contact share proportional to eligible interchain site inventory. It does not measure solvent-accessible area, affinity, or H-bond energy. Zero total contacts or zero eligible pairs leave the enrichment undefined.",
                         "Only urethane NH donors are included. Distinguish PCL and LDI ester carbonyls, urethane/urea carbonyls, and PEG ether oxygens. One-frame retention is conditional on the saved stride; it is not a residence time. Noncovalent packing and preparation history require independent seeds."]
    else:
        source=v["ldi_ftir_source"].strip()
        if not source or (v["ldi_ftirs"].strip()!=DEMO_SERIES and source==SERIES_NOTE):
            raise ValueError("Identify the origin of modified measurements and simulations.")
        r=contact_folds(v["ldi_ftirs"])
        for key,value in zip(("Recipe","Recipe + ether enrichment"),r["mae_pp"]):
            out.append(_metric(key,key+" MAE",value,"pp","Equal weight per withheld formulation."))
        result["tables"]=[_table("FTIR fitted area fraction, one formulation withheld",["Formulation","Specimens","Observed","Recipe prediction","Recipe + MD prediction"],r["rows"])]
        result["notes"]=[source,"Matching protocol: "+r["protocol"],"Target is the fitted bonded carbonyl band area fraction in a defined dry-film ATR spectrum. Record baseline, component shapes and PCL/LDI ester overlap. Independent DSC/WAXS can flag crystallinity changes. NH contact counts and FTIR band areas are different observables.",
                         "Recipe predictors are PEG mass fraction among soft segments and feed NCO/OH. Additional predictor is interchain NH–PEG ether contact enrichment. Ridge uses fixed alpha=1, standardized inputs in each training fold and an unpenalized intercept. The fitted values are unconstrained. Average replicates within each formulation before scoring. The supplied constructed series exercises the estimator; empirical performance requires matched trajectories and spectra."]
    return result
