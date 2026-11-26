import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from trajectory_analysis import analyze_replicas, box_lengths, file_hash, nearest_self_image


class ObservablesTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.d = Path(self.tmp.name)
        self.top = self.d / "topology.json"
        self.top.write_text(json.dumps(dict(
            elements=["C","O","O","H","H","O","H","H"],
            mass_da=[12,16,16,1,1,16,1,1],
            bonds=[[0,1,2],[2,3],[2,4],[5,6],[5,7]],
            solute_indices=[0,1],carbonyl_oxygens=[1], water_oxygens=[2,5])))
        self.npz = self.d / "seed_7.npz"

    def make(self, far=0.68, path=None):
        path = path or self.npz
        frames = []
        for L in [2,4,2,4]:
            # Mass-unequal dimer crosses x boundary, with both waters at known O–O separations.
            rows = np.array([[L-.1,.2,.2],[.1,.2,.2],
                             [.32,.2,.2],[.32,.27,.2],[.35,.15,.2],
                             [.1+far,.2,.2],[.13+far,.27,.2],[.1+far,.15,.2]])
            frames.append(rows)
        np.savez(path, positions_nm=frames, box_nm=np.array([np.eye(3)*L for L in [2,4,2,4]]),
                 time_ps=np.linspace(.2,.8,4))
        return path

    def test_wrapped_mass_weighted_and_snapshot(self):
        s = analyze_replicas(self.top,[self.make()],self.d/"out",bin_width_nm=.1)
        expected_rg=.2*np.sqrt(12*16)/(12+16)
        np.testing.assert_allclose(s["replicas"][0]["rg_nm"], expected_rg, atol=1e-13)
        self.assertAlmostEqual(s["rg"]["mean_nm"],expected_rg)
        self.assertIsNone(s["rg"]["between_replica_sd_nm"])
        self.assertEqual(s["coordination"]["mean_neighbors"],1)
        self.assertAlmostEqual(s["replicas"][0]["min_solute_image_heavy_distance_nm"],1.8)
        snapshot=json.loads((self.d/"out/snapshot.json").read_text())
        self.assertEqual(snapshot["selected_water_count"],1)
        atoms=snapshot["atoms"]
        self.assertAlmostEqual(atoms[1]["x"]-atoms[0]["x"],.2)
        self.assertAlmostEqual(atoms[2]["x"]-atoms[1]["x"],.22)
        self.assertEqual(snapshot["bonds"][0]["order"],2)

    def test_npt_rdf_exact_shell_and_frame_normalization(self):
        s = analyze_replicas(self.top,[self.make()],self.d/"out",bin_width_nm=.1)
        shell = (4*np.pi/3)*(.3**3-.2**3)
        self.assertAlmostEqual(s["rdf"]["mean_g"][2], np.mean([8,64,8,64])/(2*shell))
        self.assertEqual(s["rdf"]["mean_coordination"][2],1)
        self.assertEqual(s["rdf"]["mean_coordination"][-1],2)
        self.assertAlmostEqual(s["rdf"]["r_max_nm"],1)
        self.assertAlmostEqual(np.dot(np.diff(s["rg"]["bin_edges_nm"]),s["rg"]["density_per_nm"]),1)

    def test_between_replicas_not_frames_and_hydration_fixed_cutoff(self):
        other=self.make(far=.30,path=self.d/"seed_9.npz")
        s=analyze_replicas(self.top,[other,self.make()],self.d/"out",bin_width_nm=.1)
        self.assertEqual([r["seed"] for r in s["replicas"]],[7,9])
        self.assertAlmostEqual(s["coordination"]["mean_neighbors"],1.5)
        self.assertAlmostEqual(s["coordination"]["between_replica_sd"],np.sqrt(.5))
        self.assertEqual(len(s["provenance"]["outputs_sha256"]),5)
        self.assertEqual(len(s["provenance"]["trajectories_sha256"]),2)

    def test_directional_hbond_not_same_as_close_oxygen(self):
        src=self.make()
        with np.load(src) as f:
            data={k:f[k] for k in f.files}
        # Distance from donor O to acceptor O is 0.22 nm. First water H points
        # directly toward acceptor, then away; O positions never change.
        for t in range(4):
            data["positions_nm"][t,3]=[.225 if t<2 else .415,.2,.2]
        np.savez(src,**data)
        result=analyze_replicas(self.top,[src],self.d/"out")
        self.assertEqual(result["replicas"][0]["hydrogen_bonds_per_carbonyl"], [1.,1.,0.,0.])
        self.assertEqual(result["hydrogen_bonds"]["mean_per_carbonyl"],.5)
        self.assertEqual(result["coordination"]["mean_neighbors"],1.)

    def test_no_carbonyl_keeps_chain_result(self):
        top=json.loads(self.top.read_text())
        top["carbonyl_oxygens"]=[]
        self.top.write_text(json.dumps(top))
        s=analyze_replicas(self.top,[self.make()],self.d/"out",bin_width_nm=.1)
        self.assertEqual(s["rdf"]["status"],"no_matching_solute_sites")
        self.assertEqual(s["rdf"]["mean_g"],[])
        self.assertIsNone(s["coordination"]["mean_neighbors"])
        self.assertGreater(s["rg"]["mean_nm"],0)

    def test_distinct_amide_sites_reconcile_for_each_seed_and_pooled(self):
        elements = list("OCNHCCNCOCOHHOHH")
        bonds = [[0,1],[1,2],[2,3],[2,4],[4,5],[5,6],[6,7],[7,8],[6,9],
                 [10,11],[10,12],[13,14],[13,15]]
        self.top.write_text(json.dumps(dict(elements=elements,mass_da=[{"O":16,"C":12,"H":1,"N":14}[e] for e in elements],
             bonds=bonds,solute_indices=list(range(10)),carbonyl_oxygens=[0,8],water_oxygens=[10,13])))
        xyz = np.array([[.80,.5,.5],[.9,.5,.5],[1.,.5,.5],[1.,.45,.5],[1.1,.5,.5],
                        [1.25,.5,.5],[1.4,.5,.5],[1.55,.5,.5],[1.7,.5,.5],[1.4,.61,.5],
                        [.55,.5,.5],[.65,.5,.5],[.53,.59,.5],[.80,.76,.5],[.80,.66,.5],[.89,.78,.5]])
        paths = [self.npz,self.d/"seed_9.npz"]
        for seed_index,path in enumerate(paths):
            coords = xyz.copy()
            coords[13:,0] += .9*seed_index
            np.savez(path,positions_nm=np.repeat(coords[None],4,axis=0),box_nm=np.full((4,3),3.),
                     time_ps=np.arange(1,5)*.2)
        result=analyze_replicas(self.top,paths,self.d/"out")
        a,b=result["site_hydration"]
        self.assertEqual([a["amide_class"],b["amide_class"]], ["secondary_amide","tertiary_amide"])
        self.assertEqual([a["n_bonded_hydrogens"],b["n_bonded_hydrogens"]],[1,0])
        self.assertEqual(b["n_bonds_to_nearest_secondary_amide_n"],3)
        self.assertEqual([a["mean_neighbors"],b["mean_neighbors"]],[1.5,.5])
        self.assertEqual([a["mean_hbonds"],b["mean_hbonds"]],[1.5,.5])
        self.assertAlmostEqual(a["between_replica_sd_hbonds"],np.sqrt(.5))
        self.assertAlmostEqual((a["mean_neighbors"]+b["mean_neighbors"])/2,result["coordination"]["mean_neighbors"])
        self.assertAlmostEqual((a["mean_hbonds"]+b["mean_hbonds"])/2,result["hydrogen_bonds"]["mean_per_carbonyl"])
        for i,run in enumerate(result["replicas"]):
            for metric,all_sites in [("mean_hbonds","hydrogen_bond_mean_per_carbonyl"),
                                     ("mean_neighbors","coordination_mean_at_035_nm")]:
                self.assertAlmostEqual(run[all_sites],(a["replicas"][i][metric]+b["replicas"][i][metric])/2)

    def test_reject_mixed_simulation_protocol_or_missing_manifest(self):
        paths = [self.make(), self.make(path=self.d/"seed_9.npz")]
        def record(path, seed, temperature):
            path.with_suffix(".json").write_text(json.dumps({
                "parameters":{"seed":seed,"temperature_K":temperature},
                "run_id":f"run-{seed}","trajectory":path.name,"trajectory_sha256":file_hash(path),
                "shared_topology_sha256":file_hash(self.top)}))
        record(paths[0],7,298)
        with self.assertRaises(ValueError):
            analyze_replicas(self.top,paths,self.d/"out")
        record(paths[1],9,310)
        with self.assertRaises(ValueError):
            analyze_replicas(self.top,paths,self.d/"out")
        record(paths[1],9,298)
        result=analyze_replicas(self.top,paths,self.d/"out")
        self.assertEqual(result["provenance"]["simulation_run_ids"]["seed_9.npz"],"run-9")
        record(paths[1],8,298)
        with self.assertRaises(ValueError):
            analyze_replicas(self.top,paths,self.d/"out")

    def test_minimum_nonzero_image_matches_brute_force(self):
        L=np.array([2.,3.,2.5])
        xyz=np.array([[0.,0.,0.],[1.82,1.80,1.10],[-.15,.20,-.10]])
        brute=min(np.linalg.norm(a-b+L*np.asarray([i,j,k]))
                  for a in xyz for b in xyz for i in range(-2,3) for j in range(-2,3)
                  for k in range(-2,3) if (i,j,k)!=(0,0,0))
        self.assertAlmostEqual(brute,nearest_self_image(xyz,L))

    def test_reject_invalid_topology_boxes_and_repeated_seeds(self):
        triangle=np.eye(3);triangle[1,0]=.1
        with self.assertRaises(ValueError):box_lengths(np.array([triangle]))
        with self.assertRaises(ValueError):analyze_replicas(self.top,[self.make(),self.npz],self.d/"out")
        top=json.loads(self.top.read_text())
        top["bonds"]=top["bonds"][1:]
        self.top.write_text(json.dumps(top))
        with self.assertRaises(ValueError):analyze_replicas(self.top,[self.make()],self.d/"out")


if __name__ == "__main__":
    unittest.main()
