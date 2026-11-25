"""Small numerical and chemical checks for trajectory-based comparisons."""
import unittest
import json
from pathlib import Path
import numpy as np
from compare_solutions import acf_diagnostic, rank_folded_split_rhat, water_to_acceptors, solute_donor_to_water, mass_rg, retention, load_case

class TestComparison(unittest.TestCase):
    def test_directional_hydration_and_periodic_wrap(self):
        # O at 0.95, H toward 0.0, acceptor across x boundary at 0.23.
        xyz=np.array([[.95,.5,.5],[.05,.5,.5],[.95,.6,.5],[.23,.5,.5],[.67,.5,.5]])
        c,h=water_to_acceptors(xyz,np.ones(3),np.array([0]),np.array([[1,2]]),np.array([3,4]))
        np.testing.assert_array_equal(c,[1,1]);np.testing.assert_array_equal(h,[1,0])
        self.assertEqual(solute_donor_to_water(xyz,np.ones(3),np.array([3,4]),0,1),1)

    def test_leaving_and_returning_water_is_detected_only_at_endpoints(self):
        # One site, two waters; water 0 leaves and returns during lag=2,
        # while water 1 stays present through all three sampled frames.
        o=np.array([[[1,1]],[[0,1]],[[1,1]]],dtype=bool)
        den,endpoint,continuous=retention(o,2)
        np.testing.assert_array_equal(den,[2])
        np.testing.assert_array_equal(endpoint,[2])
        np.testing.assert_array_equal(continuous,[1])

    def test_original_measured_pilot_matches_published_observables(self):
        parent=Path(__file__).resolve().parents[2]
        locations=(parent/"pilot", parent/"mix08"/"pilot")
        pilot=next((p for p in locations if (p/"run/seed_13031.npz").is_file()), None)
        if pilot is None:
            self.skipTest("Source study trajectory files unavailable in standalone checkout.")
        published=json.loads((pilot/"analysis/observables.json").read_bytes())
        model,runs=load_case(pilot/"run","HH")
        self.assertEqual(model.selection()["carbonyl_oxygens"],{"alpha":3,"middle":18,"omega":14})
        for run,ref in zip(runs,published["replicas"],strict=True):
            self.assertEqual(run["seed"],ref["seed"])
            for measured,key in [("rg_nm","rg_mean_nm"),("carbonyl_neighbors","coordination_mean_at_035_nm"),
                                 ("water_to_carbonyl_hb","hydrogen_bond_mean_per_carbonyl")]:
                self.assertAlmostEqual(np.mean([f[measured] for f in run["frames"]]),ref[key],places=12)

    def test_weighted_rg_under_translation(self):
        p=np.array([[0.,0.,0.],[2.,0.,0.]])
        self.assertAlmostEqual(mass_rg(p+5,np.array([1.,3.])),np.sqrt(.75))

    def test_effective_count_and_flat_series(self):
        rng=np.random.default_rng(2164)
        trial=rng.normal(size=1000)
        independent=acf_diagnostic(trial,1.)
        correlated=acf_diagnostic(np.repeat(trial[:100],10),1.)
        self.assertGreater(independent["ess"],5*correlated["ess"])
        self.assertEqual(acf_diagnostic(np.ones(100),1.)["ess"],None)

    def test_between_chain_and_within_chain_shift(self):
        rng=np.random.default_rng(30521);x=rng.normal(size=(3,500))
        self.assertLess(rank_folded_split_rhat(x),1.03)
        x[0]+=3
        self.assertGreater(rank_folded_split_rhat(x),1.3)
        x=np.tile(np.arange(25.),(3,1))
        self.assertGreater(rank_folded_split_rhat(x),1.4)

if __name__=="__main__": unittest.main()
