from pathlib import Path
import unittest
import numpy as np
from discovery.core import RunContext
from discovery.engines.coset_group_search import coset_generators,search_coset_orbits
from discovery.engines.group_seed_search import _closure
from discovery.runner import registry


def swap(d,i,j):
 p=np.eye(d,dtype=int);p[[i,j]]=p[[j,i]]
 return {'matrix':p.tolist(),'denominator':1}

class CosetGroupTests(unittest.TestCase):
 def setUp(self):
  self.ctx=RunContext(Path('.'));self.g=[swap(4,0,1),swap(4,1,2),swap(4,2,3)];self.h=self.g[:2]
 def test_square_index_builds_exact_standard_representation(self):
  mats,info=coset_generators(self.g,self.h,4,1,0,self.ctx)
  self.assertEqual(info,{'group_order':24,'subgroup_order':6,'index':4,'standard_dimension':3,'representation_dimension':3})
  group,den=_closure(mats,3,self.ctx);self.assertEqual(len(group),24)
  for g in group:self.assertTrue(np.array_equal(g@g.T,den*den*np.eye(3,dtype=np.int64)))
 def test_direct_sum_and_trivial_coordinates(self):
  mats,info=coset_generators(self.g,self.h,4,2,1,self.ctx);self.assertEqual(info['representation_dimension'],7)
  group,den=_closure(mats,7,self.ctx);self.assertEqual(len(group),24)
  self.assertTrue(all(g[-1,-1]==den for g in group))
 def test_non_subgroup_and_nonsquare_index_rejected(self):
  with self.assertRaisesRegex(ValueError,'not contained'):coset_generators([self.g[0]],[self.g[2]],4,1,0,self.ctx)
  with self.assertRaisesRegex(ValueError,'square index'):coset_generators(self.g,self.g[:1],4,1,0,self.ctx)
 def test_supplemental_representation_checks_group_relations(self):
  mats,info=coset_generators(self.g,self.h,4,1,0,self.ctx,self.g,4);self.assertEqual(info['representation_dimension'],7)
  group,den=_closure(mats,7,self.ctx);self.assertEqual(len(group),24)
  negate={'matrix':[[-1,0],[0,-1]],'denominator':1};identity={'matrix':[[1,0],[0,1]],'denominator':1}
  with self.assertRaisesRegex(ValueError,'violate source-group relations'):
   coset_generators(self.g,self.h,4,1,0,self.ctx,[negate,identity,identity],2)
 def test_operator_registered_and_generates_provisional_search(self):
  ops,_=registry();self.assertIn('search.coset_group_orbits',ops)
  p={'dimension':4,'generators':self.g,'subgroup_generators':self.h,'copies':1,'trivial_dimensions':1,'orbit_count':1,'attempts':1,'max_iterations':5,'seconds':3,'seed':7}
  a=search_coset_orbits([],p,self.ctx);self.assertEqual(a.verification.status,'unverified');self.assertEqual(a.data['result']['coset_representation']['representation_dimension'],4)
if __name__=='__main__':unittest.main()
