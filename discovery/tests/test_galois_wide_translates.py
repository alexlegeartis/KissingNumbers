from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.galois_kerdock import construct
from discovery.engines.galois_wide_translates import search
from discovery.runner import registry

class GaloisWideTranslateTests(unittest.TestCase):
 def setUp(self):self.ctx=RunContext(Path('.'))
 def test_random_translate_search_replays_exactly(self):
  a=construct([],{'m':3},self.ctx);b=search([a],{'trials':200,'seed':7,'seconds':3,'max_layers':8,'add_d_roots':True},self.ctx)
  self.assertGreaterEqual(b.facts['translate_layers'],2);self.assertEqual(b.facts['cardinality'],256*b.facts['translate_layers']+480);self.assertEqual(b.verification.status,'exact');self.assertEqual(verify_artifact(b,self.ctx).status,'exact')
  b.data['shifts'][1]=b.data['shifts'][0][:]
  with self.assertRaisesRegex(ValueError,'invalid binary shifts'):verify_artifact(b,self.ctx)
 def test_rejects_out_of_range_shift(self):
  a=construct([],{'m':3},self.ctx);b=search([a],{'trials':1,'seed':7,'seconds':3},self.ctx);b.data['shifts']=[[1<<16,0]]
  with self.assertRaisesRegex(ValueError,'invalid binary shifts'):verify_artifact(b,self.ctx)
 def test_work_budget_returns_checked_incumbent(self):
  ctx=RunContext(Path('.'),{'max_rank_work':132000});a=construct([],{'m':3},ctx);b=search([a],{'trials':200,'seed':7,'seconds':3,'max_layers':8},ctx)
  self.assertGreaterEqual(b.facts['translate_layers'],1);self.assertEqual(verify_artifact(b,ctx).status,'exact')
 def test_registry(self):self.assertIn('galois.random_translates',registry()[0])
if __name__=='__main__':unittest.main()
