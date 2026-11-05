from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.binary_span_cosets import build
class BinarySpanCosetsTests(unittest.TestCase):
 def test_exact_union_and_duplicate_rejection(self):
  ctx=RunContext(Path('.'));a=build([] ,{'length':4,'generators':[15],'representatives':[0,3,5,6]},ctx)
  self.assertEqual(a.facts,{'length':4,'cardinality':8,'min_distance':2});self.assertEqual(verify_artifact(a,ctx).status,'exact')
  with self.assertRaises(ValueError):build([] ,{'length':4,'generators':[15],'representatives':[0,15]},ctx)
 def test_dependent_rejected(self):
  with self.assertRaises(ValueError):build([] ,{'length':4,'generators':[3,3],'representatives':[0]},RunContext(Path('.')))
if __name__=='__main__':unittest.main()
