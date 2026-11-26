"""Compute or summarize a short, seeded TIP3P NPT water reference with OpenMM."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import time

SEEDS = (13031, 27011, 49033)
N_WATER = 216
T_K, P_BAR = 298.15, 1.0
DT_PS, EQUIL_PS, PROD_PS, SAMPLE_PS = 0.002, 12, 32, 0.2
BLOCK_PS = 8
BASE = Path(__file__).resolve().parent

def digest(data):
    return hashlib.sha256(data).hexdigest()

def json_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()

def run(seed, out):
    import openmm as mm
    from openmm import app, unit
    from importlib.metadata import version
    import numpy as np

    start = time.monotonic()
    out.mkdir(parents=True, exist_ok=True)
    forcefield = app.ForceField("tip3p.xml")
    model = app.Modeller(app.Topology(), [])
    model.addSolvent(forcefield, model="tip3p", numAdded=N_WATER, neutralize=False)
    system = forcefield.createSystem(model.topology, nonbondedMethod=app.PME,
                                     nonbondedCutoff=0.8*unit.nanometer, rigidWater=True,
                                     ewaldErrorTolerance=0.0005)
    barostat = mm.MonteCarloBarostat(P_BAR*unit.bar, T_K*unit.kelvin, 25)
    barostat.setRandomNumberSeed(seed+2)
    system.addForce(barostat)
    integrator = mm.LangevinMiddleIntegrator(T_K*unit.kelvin, 1/unit.picosecond,
                                             DT_PS*unit.picoseconds)
    integrator.setRandomNumberSeed(seed+1)
    integrator.setConstraintTolerance(1e-5)
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    sim = app.Simulation(model.topology, system, integrator, mm.Platform.getPlatformByName("CPU"),
                         {"Threads": "1", "DeterministicForces": "true"})
    sim.context.setPositions(model.positions)
    sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer, maxIterations=200)
    sim.context.setVelocitiesToTemperature(T_K*unit.kelvin, seed)
    volume_conversion = 602.214076  # nm3 × g/mL × N_A -> daltons
    mass_da = sum(system.getParticleMass(i).value_in_unit(unit.dalton)
                  for i in range(system.getNumParticles()))
    xml = (Path(app.__file__).parent/"data"/"tip3p.xml").read_bytes()
    protocol = {
        "model": "Single-component periodic cubic box; 216 rigid neutral TIP3P waters",
        "engine": "OpenMM", "engine_version": version("openmm"), "engine_build": mm.version.version,
        "platform": "CPU", "threads_per_replica": 1, "deterministic_forces": True,
        "temperature_K": T_K, "pressure_bar": P_BAR, "ensemble": "isothermal-isobaric (NPT)",
        "force_field": "openmm.app/data/tip3p.xml", "force_field_sha256": digest(xml),
        "n_water": N_WATER, "n_atoms": system.getNumParticles(), "n_constraints": system.getNumConstraints(),
        "water_model": "TIP3P", "component_smiles": "O", "ionic_strength_mol_l": 0.0,
        "ensemble_formula": "mass_da / (602.214076 * volume_nm3)",
        "integrator": "LangevinMiddleIntegrator", "friction_per_ps": 1, "step_ps": DT_PS,
        "barostat": "isotropic MonteCarloBarostat", "barostat_frequency_steps": 25,
        "pme": True, "nonbonded_cutoff_nm": 0.8, "ewald_relative_tolerance": 0.0005,
        "lj_dispersion_correction": nb.getUseDispersionCorrection(), "lj_switch": nb.getUseSwitchingFunction(),
        "remove_com_frequency_steps": next(f.getFrequency() for f in system.getForces() if isinstance(f, mm.CMMotionRemover)),
        "constraint_tolerance": integrator.getConstraintTolerance(),
        "minimization_tolerance_kj_mol_nm": 10, "minimization_max_iterations": 200,
        "common_start": "OpenMM Modeller.addSolvent(numAdded=216, model='tip3p')",
        "starting_positions_sha256": digest(np.asarray(model.positions.value_in_unit(unit.nanometer),dtype="<f8").tobytes()),
        "initial_cube_edge_nm": model.topology.getUnitCellDimensions().value_in_unit(unit.nanometer)[0],
        "equilibration_ps": EQUIL_PS, "production_ps": PROD_PS, "sample_interval_ps": SAMPLE_PS,
        "production_block_size_ps": BLOCK_PS, "mass_da": mass_da, "seeds": SEEDS,
        "seed_roles": "initial velocity=seed; Langevin=seed+1; barostat=seed+2",
    }
    dof = system.getNumParticles()*3 - system.getNumConstraints()-3
    rows = []
    sample_steps = round(SAMPLE_PS/DT_PS)
    for stage, length in (("equilibration",EQUIL_PS), ("production",PROD_PS)):
        for i in range(1,round(length/SAMPLE_PS)+1):
            sim.step(sample_steps)
            state=sim.context.getState(energy=True)
            vol=state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
            pot=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            kinetic=state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            rho=mass_da/(volume_conversion*vol)
            if not all(math.isfinite(x) for x in (vol,pot,kinetic,rho)) or vol<=0:
                raise ValueError("Nonfinite NPT sample")
            rows.append({"seed":seed,"stage":stage,"stage_time_ps":round(i*SAMPLE_PS,5),
                         "volume_nm3":round(vol,7),"density_g_ml":round(rho,8),
                         "cell_potential_energy_kj_mol":round(pot,5),
                         "kinetic_temperature_K":round(2*kinetic/(dof*0.00831446261815324),5)})
        print(seed,stage,"complete",round(time.monotonic()-start,1),"s",flush=True)
    # State carries full-precision periodic box, positions and velocities for geometric reuse.
    # A new solute or force field requires its own preparation, sampling and equilibration.
    final=sim.context.getState(positions=True, velocities=True, enforcePeriodicBox=True, parameters=True)
    state_name=f"water_{seed}_state.xml.gz"
    (out/state_name).write_bytes(gzip.compress(mm.XmlSerializer.serialize(final).encode(),mtime=0))
    dump={"seed":seed,"protocol":protocol,"rows":rows,"final_state":state_name,
          "state_sha256":digest((out/state_name).read_bytes())}
    (out/f"water_{seed}_run.json").write_bytes(json_bytes(dump))

def summarize(out):
    runs=[json.loads((out/f"water_{s}_run.json").read_bytes()) for s in SEEDS]
    protocol=runs[0]["protocol"]
    if any(r["protocol"] != protocol for r in runs):
        raise ValueError("Replica protocols or common starting structures differ")
    model_id="sha256:"+digest(json_bytes({"schema":"water-npt/1","protocol":protocol}))
    reps=[]
    with (out/"water_npt.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=("seed","stage","stage_time_ps","volume_nm3","density_g_ml","cell_potential_energy_kj_mol","potential_energy_per_water_kj_mol","kinetic_temperature_K"))
        writer.writeheader()
        for r in runs:
            for row in r["rows"]:
                row["potential_energy_per_water_kj_mol"] = round(row["cell_potential_energy_kj_mol"]/protocol["n_water"],8)
            writer.writerows(r["rows"])
            prod=[p for p in r["rows"] if p["stage"]=="production"]
            if len(prod)!=round(PROD_PS/SAMPLE_PS):
                raise ValueError("Incomplete NPT trajectory")
            rho=[p["density_g_ml"] for p in prod]
            step=round(BLOCK_PS/SAMPLE_PS)
            blocks=[round(statistics.mean(rho[i:i+step]),8) for i in range(0,len(rho),step)]
            reps.append({"seed":r["seed"],"equilibration_ps":EQUIL_PS,"production_ps":PROD_PS,
                         "density_mean_g_ml":round(statistics.mean(rho),8),"block_means_g_ml":blocks,
                         "potential_energy_per_water_kj_mol":round(statistics.mean(p["potential_energy_per_water_kj_mol"] for p in prod),8),
                         "density_timeseries":[{"time_ps":p["stage_time_ps"],"density_g_ml":p["density_g_ml"]} for p in prod],
                         "final_state":r["final_state"],"final_state_sha256":r["state_sha256"]})
    means=[r["density_mean_g_ml"] for r in reps]
    energies=[r["potential_energy_per_water_kj_mol"] for r in reps]
    report={"generator_sha256":digest(Path(__file__).read_bytes()),"schema":"water-npt/1","id":model_id,"phase":"pure_liquid","component_smiles":"O",
            "engine":protocol["engine"],"engine_version":protocol["engine_version"],"water_model":"TIP3P",
            "temperature_K":T_K,"pressure_bar":P_BAR,"replicates":reps,
            "density_mean_g_ml":round(statistics.mean(means),8),
            "potential_energy_per_water_kj_mol":round(statistics.mean(energies),8),
            "replica_mean_potential_sd_kj_per_mol_water":round(statistics.stdev(energies),8),
            "potential_energy_method":"Mean classical cell potential divided by 216 waters over production samples; kJ per mol H2O. Specific to this potential, cutoff, electrostatic conventions and state. No free energy or heat of mixing is inferred.",
            "replica_mean_sd_g_ml":round(statistics.stdev(means),8),
            "quality_status":"preliminary","protocol":protocol,
            "source_script":"run_water.py","trajectory_csv":"water_npt.csv","trajectory_sha256":digest((out/"water_npt.csv").read_bytes()),
            "uncertainty":"Sample standard deviation across three replica means; descriptive scatter. Four 8 ps blocks per replica; no independence or convergence guarantee.",
            "reuse":"Final periodic geometries for pure TIP3P water; record model and protocol. Solute insertion changes composition and requires new parameterization and equilibration."}
    (out/"water_npt.json").write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("Density",report["density_mean_g_ml"],"g/mL; SD of replica means",report["replica_mean_sd_g_ml"])

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task",choices=("run","summarize","all"))
    parser.add_argument("--seed",type=int,choices=SEEDS)
    parser.add_argument("--out",type=Path,default=BASE)
    args=parser.parse_args()
    if args.task=="run" and args.seed is None:
        parser.error("run requires --seed")
    os.environ.setdefault("OPENMM_CPU_THREADS","1")
    if args.task=="all":
        for seed in SEEDS: run(seed,args.out)
    elif args.task=="run": run(args.seed,args.out)
    if args.task in ("all","summarize"): summarize(args.out)
