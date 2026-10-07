import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from targetscan.collision_floor import collision_mse_floor as floor

class FloorTests(unittest.TestCase):
    def test_duplicate_pair(self):
        self.assertEqual(floor([b'a',b'a'],[1,3])['empirical_mse_floor'],1.)
    def test_unique_zero(self):
        self.assertEqual(floor([b'a',b'b'],[1,3])['empirical_mse_floor'],0.)
    def test_weighted_row_count(self):
        self.assertAlmostEqual(floor([b'a',b'a',b'b'],[1,3,99])['empirical_mse_floor'],2/3)
    def test_common_translation(self):
        self.assertEqual(floor([b'a',b'a'],[100001,100003])['empirical_mse_floor'],1.)
    def test_permutation(self):
        self.assertEqual(floor([b'a',b'a',b'b'],[1,3,99]),floor([b'b',b'a',b'a'],[99,3,1]))
    def test_reject_empty_and_length(self):
        for k,y in [([],[]),([b'a'],[])]:
            with self.assertRaises(ValueError):floor(k,y)
    def test_reject_nonfinite(self):
        for v in [float('inf'),float('nan')]:
            with self.assertRaises(ValueError):floor([b'a'],[v])
    def test_reject_coercion(self):
        for v in [True,'3',None]:
            with self.assertRaises(TypeError):floor([b'a'],[v])
        with self.assertRaises(TypeError):floor(['id'],[3])
    def test_range_rejected(self):
        with self.assertRaises(ValueError):floor([b'a',b'a'],[-1e308,1e308])
    def test_oracle_vs_constant(self):
        for a in range(-5,6):
            for b in range(-5,6):
                f=floor([b'k',b'k'],[a,b])['empirical_mse_floor']
                for p in range(-7,8):self.assertLessEqual(f,((a-p)**2+(b-p)**2)/2+1e-12)

if __name__=='__main__':unittest.main()
