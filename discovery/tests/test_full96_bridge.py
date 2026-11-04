import copy
import unittest
from discovery.benchmarks.run_full96 import _compare_reports,EXPECTED_96

class Full96BridgeTests(unittest.TestCase):
    def fixture(self):
        terms={'dense':137438953472,'transversal24':4294967296,'support24':426344448,'support6':2097152,'pure':192}
        a={'status':'PASS','support_sha256':'synthetic-comparison-fixture','dimensions':{str(d):{} for d in range(88,97)}}
        a['dimensions']['96']={'total':EXPECTED_96,'terms':terms}
        return a
    def test_requires_full_layers_and_independent_agreement(self):
        a=self.fixture();self.assertEqual(_compare_reports(a,copy.deepcopy(a))['total'],EXPECTED_96)
        for mutate in (lambda x:x['dimensions']['96']['terms'].pop('support6'),lambda x:x['dimensions'].pop('88'),lambda x:x.update(status='FAIL')):
            bad=copy.deepcopy(a);mutate(bad)
            with self.assertRaises(ValueError):_compare_reports(bad,copy.deepcopy(bad))
        bad=copy.deepcopy(a);bad['dimensions']['96']['total']+=1
        with self.assertRaises(ValueError):_compare_reports(a,bad)
