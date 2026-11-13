"""Small exact-arithmetic tests for the 26D/27D search helpers."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

HERE=Path(__file__).parent

def load(stem):
    spec=importlib.util.spec_from_file_location(stem,HERE/(stem+'.py'))
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

class Test26D(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.m=load('solve_joint26')
    def test_qsqrt3_predicate(self):
        for a in range(-20,41):
            for b in range(-30,31):
                got=bool(self.m.above16(np.array([a]),np.array([b]))[0])
                self.assertEqual(got,a+b*np.sqrt(3)>16)
    def test_components(self):
        blockers=[np.array([7]),np.array([7]),np.array([8])]
        cs=self.m.components(3,blockers,np.array([0]),np.array([2]))
        self.assertEqual(len(cs),1)
        self.assertEqual(cs[0][0].tolist(),[0,1,2])
    def test_shared_blocker_greedy(self):
        blocks=[np.array([7]),np.array([7]),np.array([8])]
        gain,chosen,deleted=self.m.greedy(np.arange(3),blocks,{(0,2)},7)
        self.assertEqual(gain,1)
        self.assertEqual(set(deleted),{7})

class Test27D(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.m=load('search_partial_fibers27')
    def test_catalogues(self):
        H=self.m.heads2160(); T=self.m.tails12()
        self.assertEqual(H.shape,(2160,8)); self.assertEqual(T.shape,(12,3))
        self.assertTrue(np.all((H.astype(np.int64)**2).sum(1)==16))
        self.assertTrue(np.all((T.astype(np.int64)**2).sum(1)==8))
    def test_qsqrt2_predicate(self):
        for p in range(-8,33):
            for q in range(-12,13):
                got=bool(self.m.greater2(np.array([p]),np.array([q]))[0])
                self.assertEqual(got,p+np.sqrt(2)*q>16)

if __name__=='__main__': unittest.main()
