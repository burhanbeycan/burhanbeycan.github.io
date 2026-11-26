import unittest
import numpy as np

from amide_torsions import torsions, verify


class DihedralDefinition(unittest.TestCase):
    def test_signed_quarter_turns_across_periodic_boundary(self):
        # The middle bond is along +x. Reversing the last bond's z orientation
        # should reverse the angle, including when atoms are wrapped into [0,L).
        x=np.array([[0,1,0],[0,0,0],[1,0,0],[1,0,1]],float)
        p=np.stack([x,x*np.array([1,1,-1])])
        wrapped=(p+np.array([3.6,3.5,3.8]))%4
        b=np.tile(np.eye(3)[None,:,:]*4,(2,1,1))
        a,c=torsions(wrapped,b,[0,1,2,3])
        np.testing.assert_allclose(a,[90,-90],atol=1e-10)
        verify(c,a)


if __name__=='__main__':
    unittest.main()
