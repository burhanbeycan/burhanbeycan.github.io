"""Verify stored MD output without an MD runtime or a trajectory rerun."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import statistics
import unittest
from xml.etree import ElementTree

BASE = Path(__file__).resolve().parent

class StoredWaterRecord(unittest.TestCase):
    def test_model_and_replicates_match_raw_frames(self):
        report = json.loads((BASE / "water_npt.json").read_text())
        protocol = report["protocol"]
        encoded = json.dumps({"schema": report["schema"], "protocol": protocol},
                             sort_keys=True, separators=(",", ":")).encode()
        self.assertEqual(report["id"], "sha256:" + hashlib.sha256(encoded).hexdigest())
        self.assertEqual(report["generator_sha256"], hashlib.sha256((BASE / "run_water.py").read_bytes()).hexdigest())
        csv_path = BASE / report["trajectory_csv"]
        self.assertEqual(report["trajectory_sha256"], hashlib.sha256(csv_path.read_bytes()).hexdigest())
        with csv_path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 660)
        means, energies = [], []
        for replica in report["replicates"]:
            seed = replica["seed"]
            prod = [r for r in rows if r["stage"]=="production" and int(r["seed"])==seed]
            self.assertEqual(len(prod), 160)
            rho = [float(row["density_g_ml"]) for row in prod]
            for row, point in zip(prod,replica["density_timeseries"],strict=True):
                self.assertEqual(point, {"time_ps":float(row["stage_time_ps"]), "density_g_ml":float(row["density_g_ml"])})
                self.assertAlmostEqual(point["density_g_ml"],protocol["mass_da"]/(602.214076*float(row["volume_nm3"])),places=6)
                self.assertAlmostEqual(float(row["potential_energy_per_water_kj_mol"]),float(row["cell_potential_energy_kj_mol"])/protocol["n_water"],places=7)
            means.append(statistics.mean(rho))
            energies.append(statistics.mean(float(r["potential_energy_per_water_kj_mol"]) for r in prod))
            self.assertAlmostEqual(replica["density_mean_g_ml"],means[-1],places=7)
            self.assertAlmostEqual(replica["potential_energy_per_water_kj_mol"],energies[-1],places=7)
            self.assertEqual(replica["block_means_g_ml"],[round(statistics.mean(rho[i:i+40]),8) for i in range(0,160,40)])
            raw = (BASE / replica["final_state"]).read_bytes()
            self.assertEqual(hashlib.sha256(raw).hexdigest(),replica["final_state_sha256"])
            state = ElementTree.fromstring(gzip.decompress(raw))
            self.assertEqual(len(state.find("Positions")),protocol["n_atoms"])
            self.assertEqual(len(state.find("Velocities")),protocol["n_atoms"])
            self.assertEqual(len(state.find("PeriodicBoxVectors")),3)
        self.assertEqual(len(means),3)
        self.assertAlmostEqual(report["density_mean_g_ml"],statistics.mean(means),places=7)
        self.assertAlmostEqual(report["potential_energy_per_water_kj_mol"],statistics.mean(energies),places=7)
        self.assertAlmostEqual(report["replica_mean_sd_g_ml"],statistics.stdev(means),places=7)

if __name__ == "__main__":
    unittest.main()
