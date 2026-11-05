from pathlib import Path
import unittest
from discovery.core import Artifact, RunContext
from discovery.certification import verify_artifact
from discovery.engines.partial_orbits import rational_image_exchanges

class RationalImageTests(unittest.TestCase):
 def setUp(self):
  self.c=RunContext(Path('.'))
  self.src=Artifact('spherical_code',{'coordinates':[[1,0],[-1,0]],'weights':[1,1],'norm_squared':1})
 def call(self,gens,**kw):
  return rational_image_exchanges([self.src],{'generators':gens,'seconds':2,**kw},self.c)
 def test_registered(self):
  from discovery.runner import registry
  operators,kinds=registry();self.assertIn('search.rational_image_exchanges',operators)
  self.assertEqual(kinds['search.rational_image_exchanges'],('spherical_code',))
 def test_exchange_finds_four(self):
  s={'numerator':[[0,1],[1,0]],'denominator':1}
  result=self.call([s]);self.assertEqual(result.facts['cardinality'],4);self.assertEqual(result.data['search']['gain'],2)
  self.assertEqual(verify_artifact(result,self.c).status,'exact')
 def test_noncommuting_words(self):
  x={'numerator':[[-1,0],[0,1]],'denominator':1};s={'numerator':[[0,1],[1,0]],'denominator':1}
  result=self.call([x,s],words=[[0,1],[1,0]])
  self.assertEqual(result.facts['cardinality'],4);self.assertEqual(result.data['search']['group_words'],[[0,1],[1,0]])
 def test_word_order_is_not_reversed(self):
  one=Artifact('spherical_code',{'coordinates':[[1,0]],'weights':[1,1],'norm_squared':1})
  x={'numerator':[[-1,0],[0,1]],'denominator':1};s={'numerator':[[0,1],[1,0]],'denominator':1}
  a=rational_image_exchanges([one],{'generators':[x,s],'words':[[0,1]],'seconds':2},self.c)
  b=rational_image_exchanges([one],{'generators':[x,s],'words':[[1,0]],'seconds':2},self.c)
  self.assertEqual(set(map(tuple,a.data['coordinates'])),{(1,0),(0,-1)})
  self.assertEqual(set(map(tuple,b.data['coordinates'])),{(1,0),(0,1)})
 def test_nonintegral_rotation(self):
  rot={'numerator':[[3,-4],[4,3]],'denominator':5};result=self.call([rot])
  self.assertEqual(result.data['search']['common_denominator'],5)
  self.assertEqual(result.data['norm_squared'],25);self.assertEqual(verify_artifact(result,self.c).status,'exact')
 def test_invalid_and_bad_words(self):
  with self.assertRaises(ValueError):self.call([{'numerator':[[1,1],[0,1]],'denominator':1}])
  with self.assertRaises(ValueError):self.call([{'numerator':[[1,0],[0,1]],'denominator':1}],words=[[1]])
 def test_weighted_isometry(self):
  a=Artifact('spherical_code',{'coordinates':[[2,0],[-2,0]],'weights':[1,4],'norm_squared':4})
  swap={'numerator':[[0,4],[1,0]],'denominator':2}
  b=rational_image_exchanges([a],{'generators':[swap],'seconds':2},self.c)
  self.assertEqual(b.facts['cardinality'],4);self.assertEqual(verify_artifact(b,self.c).status,'exact')
 def test_duplicate_and_empty_words(self):
  s={'numerator':[[0,1],[1,0]],'denominator':1}
  b=self.call([s],words=[[],[0],[0]])
  self.assertEqual(b.facts['cardinality'],4);self.assertEqual(b.data['search']['distinct_new_orbits'],2)
if __name__=='__main__':unittest.main()
