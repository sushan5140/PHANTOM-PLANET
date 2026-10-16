import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from phantom.sim import simulate
from phantom.baseline import periodic_candidate
import numpy as np
class TestP01(unittest.TestCase):
    def test_deterministic(self):
        a=simulate('planet',3); b=simulate('planet',3)
        for x,y in zip(a,b): np.testing.assert_array_equal(x,y)
    def test_ground_truth(self):
        for case in ('planet','activity','both','noise'):
            t,rv,err,p,a=simulate(case)
            self.assertEqual(len(t),90)
            self.assertTrue(np.all(err>0))
            self.assertEqual(bool(np.any(p)),case in ('planet','both'))
            self.assertEqual(bool(np.any(a)),case in ('activity','both'))
    def test_periodogram(self):
        t,rv,_,_,_=simulate('planet',5)
        result=periodic_candidate(t,rv)
        self.assertGreater(result['period_days'],3)
        self.assertLess(result['period_days'],45)
    def test_bad_case(self):
        with self.assertRaises(ValueError):simulate('invalid')
if __name__=='__main__':unittest.main()
