from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.binary_grouped_cosets import search,certify
from discovery.runner import registry
class BinaryGroupedCosetsTests(unittest.TestCase):
 def test_discovery_and_independent_certification(self):
  ctx=RunContext(Path('.'));p={'length':4,'parent_generator':[3,12],'codimension':1,'subgroup_dimension':0,'minimum_distance':2,'subgroup_trials':1,'seconds':5,'max_candidates':16}
  a=search([],p,ctx);self.assertIsNotNone(a.data['proposal']);self.assertEqual(verify_artifact(a,ctx).status,'unverified');b=certify([a],{},ctx);self.assertGreater(b.facts['cardinality'],4);self.assertEqual(b.facts['min_distance'],2);self.assertEqual(verify_artifact(b,ctx).status,'exact')
 def test_registry_and_invalid_subgroup(self):
  r,_=registry();self.assertIn('search.binary_grouped_cosets',r)
  p={'length':4,'parent_generator':[3,12],'codimension':1,'subgroup_dimension':2,'minimum_distance':2,'seconds':5,'max_candidates':16}
  with self.assertRaises(ValueError):search([],p,RunContext(Path('.')))
if __name__=='__main__':unittest.main()
