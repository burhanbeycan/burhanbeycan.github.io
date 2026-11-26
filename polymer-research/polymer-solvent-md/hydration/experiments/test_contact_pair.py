"""Known-geometry controls for the two-chain simulation and contact descriptors."""
import unittest
import numpy as np
import openmm as mm

from contact_analysis import event_lists, interposed_waters, interposition_refs
from contact_pair import duplicate_solute, validate_pair, make_contact


class PairControls(unittest.TestCase):
    def test_native_pair_term_fidelity_and_no_cross_exclusions(self):
        s=mm.System()
        nb=mm.NonbondedForce();s.addForce(nb)
        bond=mm.HarmonicBondForce();s.addForce(bond)
        ang=mm.HarmonicAngleForce();s.addForce(ang)
        tor=mm.PeriodicTorsionForce();s.addForce(tor)
        for q in (-.2,.13,.12,-.05):
            s.addParticle(12.01)
            nb.addParticle(q,.32,.35)
        s.addConstraint(0,1,.145)
        bond.addBond(2,3,.151,12600)
        ang.addAngle(0,1,2,1.9,41)
        tor.addTorsion(0,1,2,3,3,.4,1.25)
        nb.addException(0,1,0,1,0);nb.addException(0,3,.008,.30,.20)
        double=duplicate_solute(s)
        self.assertEqual(validate_pair(s,double)['native_atoms'],4)
        nb2=next(f for f in double.getForces() if isinstance(f,mm.NonbondedForce))
        self.assertEqual(nb2.getNumExceptions(),4)
        for i in range(4):
            a,b,*_=nb2.getExceptionParameters(i)
            self.assertEqual(a//4,b//4)
        nb2.setParticleParameters(4,0,.32,.35)
        with self.assertRaisesRegex(ValueError,'mismatch'):
            validate_pair(s,double)

    def test_interposition_requires_original_opening_and_opposite_sites(self):
        L=np.eye(3)*2.0
        # Sites across the x boundary, initially 0.40 nm apart.
        old=np.array([[1.8,.5,.5],[.2,.5,.5],[1,.5,.5],[0,.5,.5]],dtype=float)
        refs=interposition_refs(old,L,[0],[1])
        self.assertEqual(len(refs),1)
        self.assertFalse(interposed_waters(old,L,refs,[3],0.,1))
        cur=old.copy();cur[0,0]=1.71;cur[1,0]=.29
        event=interposed_waters(cur,L,refs,[3],1.,1)[0]
        self.assertAlmostEqual(event['current_ab_nm'],.58)
        self.assertAlmostEqual(event['angle_a_ow_b_deg'],180.)
        self.assertEqual(event['water_oxygen'],3)
        # Same-site shell counts can both pass, while the angle fails.
        cur[0,0]=1.772;cur[1,0]=.228;cur[3,1]=.75
        self.assertFalse(interposed_waters(cur,L,refs,[3],1.,1))

    def test_censored_shared_water_episodes_keep_distinct_identities(self):
        times=np.arange(5.)
        flags=np.array([[1,0],[1,0],[0,1],[1,1],[0,1]],dtype=bool)
        episodes=event_lists(times,flags,[103,106])
        self.assertEqual([(v['water_oxygen_index'],v['sampled_frames']) for v in episodes],[(103,2),(106,3),(103,1)])
        self.assertTrue(episodes[0]['left_censored'])
        self.assertTrue(episodes[1]['right_censored'])


if __name__=='__main__':
    unittest.main()
