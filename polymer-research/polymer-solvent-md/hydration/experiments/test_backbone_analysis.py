import unittest
import numpy as np
from backbone_analysis import geometric_series


class TestMatchedGeometry(unittest.TestCase):
    def test_unwrap_nonuniform_mass_and_endpoint(self):
        n=9
        x=np.arange(n)*.15
        mass=np.array([14.,12.,12.]*3)
        positions=np.zeros((2,n,3)); positions[:,:,0]=x
        positions[1,:,0]=(x+1.8)%2
        boxes=np.array([np.eye(3)*2]*2)
        rg,span=geometric_series(positions,boxes,tuple(range(n)),mass)
        expected=np.sqrt(np.average((x-np.average(x,weights=mass))**2,weights=mass))
        np.testing.assert_allclose(rg,[expected]*2,atol=1.e-12)
        np.testing.assert_allclose(span,[1.2]*2)
        # End-to-end distance exceeds half the box; endpoint-only minimum-image wrapping would give 0.8 nm.
        self.assertGreater(span[1],1.)

    def test_invalid_bond_detected(self):
        coord=np.zeros((1,9,3))
        with self.assertRaisesRegex(ValueError,'Unphysical'):
            geometric_series(coord,[np.eye(3)*2],tuple(range(9)),[12.]*9)


if __name__=='__main__':unittest.main()
