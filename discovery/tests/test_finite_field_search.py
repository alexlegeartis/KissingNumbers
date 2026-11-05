from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.engines.finite_field_search import affine_flats,generate,exchange
from discovery.runner import registry

class FiniteFieldSearchTests(unittest.TestCase):
 def setUp(self):self.ctx=RunContext(Path('.'))
 def test_affine_flats_and_generation(self):
  flats=affine_flats(3,2,self.ctx);self.assertEqual(len(flats),14);self.assertEqual(len(set(flats)),14)
  p={'m':3,'modulus':11,'flat_rank':2,'terms':[[1,1]],'translations':True};r=generate([],p,self.ctx)
  self.assertEqual(len(r.data['supports']),14);self.assertTrue(all(x.bit_count()==4 for x in r.data['supports']))
  self.assertEqual(r.verification.status,'unverified')
 def test_rejects_bad_field_and_nonpermutation(self):
  p={'m':4,'modulus':21,'flat_rank':3,'terms':[[1,1]]}
  with self.assertRaisesRegex(ValueError,'reducible'):generate([],p,self.ctx)
  p={'m':4,'modulus':19,'flat_rank':3,'terms':[[1,0]]}
  with self.assertRaisesRegex(ValueError,'permutation'):generate([],p,self.ctx)
 def test_exchange_constructs_unverified_line_proposal(self):
  p={'m':4,'modulus':19,'flat_rank':3,'terms':[[1,1]]};r=generate([],p,self.ctx)
  a=exchange([r],{'seed_supports':[],'exceptional_lines':[],'parities':[0,1],'candidate_limit':1,'seconds':5},self.ctx)
  self.assertEqual(a.facts['gain_lines'],64);self.assertEqual(a.facts['gain_points_if_certified'],128);self.assertEqual(a.verification.status,'unverified')
  lines=a.data['candidate_lines'];self.assertEqual(len(lines),64)
  for i,(s,v) in enumerate(lines):
   for t,w in lines[:i]:self.assertLessEqual(abs((s&t).bit_count()-2*((s&t)&(v^w)).bit_count()),4)
 def test_overlapping_full_source_cannot_increase_single_fiber(self):
  p={'m':4,'modulus':19,'flat_rank':3,'terms':[[1,1]]};r=generate([],p,self.ctx);s=r.data['supports'][0]
  old=(s & ((1<<5)-1))
  inside=[i for i in range(16) if s>>i&1];outside=[i for i in range(16) if not(s>>i&1)]
  t=sum(1<<i for i in inside[:5]+outside[:3]);self.assertEqual((s&t).bit_count(),5)
  a=exchange([r],{'seed_supports':[t],'exceptional_lines':[],'parities':[0],'candidate_limit':1,'seconds':5},self.ctx)
  self.assertEqual(a.facts['gain_lines'],0)
 def test_coefficient_grid_skips_nonpermutations(self):
  p={'m':3,'modulus':11,'flat_rank':2,'exponents':[1,2],'coefficient_options':[[1],[0,1]]}
  r=generate([],p,self.ctx);self.assertEqual(r.data['polynomials_examined'],2);self.assertEqual(r.data['permutations'],1);self.assertEqual(r.data['nonpermutations'],1);self.assertEqual(len(r.data['supports']),14)
 def test_resource_budget_and_invalid_ambiguous_polynomial(self):
  p={'m':3,'modulus':11,'flat_rank':2,'terms':[[1,1]],'exponents':[1],'coefficient_options':[[1]]}
  with self.assertRaisesRegex(ValueError,'provide terms or'):generate([],p,self.ctx)
  p={'m':3,'modulus':11,'flat_rank':2,'exponents':[1],'coefficient_options':[[0,1,2]]}
  with self.assertRaisesRegex(ValueError,'max_trials'):generate([],p,RunContext(Path('.'),{'max_trials':2}))
 def test_sparse_coordinate_embeddings_in_96d(self):
  p={'m':3,'modulus':11,'flat_rank':3,'terms':[[1,1]],'dimension':96,'coordinate_maps':[list(range(8)),list(range(88,96))]}
  r=generate([],p,self.ctx);self.assertEqual(set(r.data['supports']),{(1<<8)-1,((1<<8)-1)<<88})
  p['coordinate_maps']=[[0]*8]
  with self.assertRaisesRegex(ValueError,'injections'):generate([],p,self.ctx)
 def test_reuses_exceptional_support_without_readding_its_line(self):
  p={'m':4,'modulus':19,'flat_rank':3,'terms':[[1,1]]};pool=generate([],p,self.ctx);s=pool.data['supports'][0]
  r=exchange([pool],{'seed_supports':[],'exceptional_lines':[[s,0]],'parities':[0],'candidate_limit':1,'seconds':5},self.ctx)
  self.assertEqual(r.facts['gain_lines'],63);self.assertNotIn([s,0],r.data['candidate_lines'])
 def test_registry(self):
  ops,_=registry();self.assertIn('finite_field.support_pool',ops);self.assertIn('finite_field.octad_exchange',ops)
if __name__=='__main__':unittest.main()
