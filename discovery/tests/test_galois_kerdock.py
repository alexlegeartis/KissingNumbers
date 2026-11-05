from pathlib import Path
import unittest
import numpy as np
from discovery.core import Artifact,RunContext
from discovery.certification import verify_artifact
from discovery.engines.galois_kerdock import _construct,construct,spherical,translate_search,search_even_extensions,find_translates
from discovery.runner import registry

class GaloisKerdockTests(unittest.TestCase):
 def setUp(self):self.ctx=RunContext(Path('.'))
 def test_nordstrom_robinson_trace_and_nonlinearity(self):
  a=construct([],{'m':3},self.ctx);self.assertEqual(a.facts['binary_length'],16);self.assertEqual(a.facts['cardinality'],256);self.assertEqual(a.facts['minimum_distance'],6);self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
  data,facts,W=_construct(3,self.ctx);gray=np.asarray(data['gray_map'],np.uint8)[W];bits=np.empty((256,16),np.uint8);bits[:,0::2]=(gray>>1)&1;bits[:,1::2]=gray&1;packed=np.packbits(bits,axis=1,bitorder='little');words={bytes(x) for x in packed};self.assertNotIn(bytes(packed[1]^packed[4]),words)
 def test_kerdock64_exact_distance(self):
  a=construct([],{'m':5},self.ctx);self.assertEqual(a.facts['binary_length'],64);self.assertEqual(a.facts['cardinality'],4096);self.assertEqual(a.facts['minimum_distance'],28);self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
 def test_spherical_embedding_and_puncture_validation(self):
  a=construct([],{'m':3},self.ctx);r=spherical([a],{},self.ctx);self.assertEqual(r.facts['cardinality'],256);self.assertEqual(verify_artifact(r,self.ctx).status,'exact')
  with self.assertRaisesRegex(ValueError,'distinct binary positions'):spherical([a],{'coordinates':[0,0]},self.ctx)
 def test_tamper_and_degree_rejected(self):
  a=construct([],{'m':3},self.ctx);a.data['generator'][0][0]=(a.data['generator'][0][0]+1)%4
  with self.assertRaisesRegex(ValueError,'modified'):verify_artifact(a,self.ctx)
  with self.assertRaisesRegex(ValueError,'degrees'):construct([],{'m':4},self.ctx)
 def test_translate_bridge(self):
  a=construct([],{'m':3},self.ctx);p={'alphabet':[[0] for _ in range(8)],'references':[[0]*8],'minimum_distance':4,'rounds':1,'oracle_steps':2,'exact_pair_limit':65536,'seconds':3};r=translate_search([a],p,self.ctx);self.assertEqual(r.data['backend'],'binary_translate_search');self.assertEqual(r.verification.status,'unverified');self.assertEqual(r.data['result']['parameters']['generator'],a.data['generator'])
 def test_searches_96d_nonlinear_evaluation_extension(self):
  a=construct([],{'m':5},self.ctx);b=search_even_extensions([a],{'count':16,'trials':100,'seed':15,'seconds':10},RunContext(Path('.'),{'solver_seconds':15}));self.assertEqual(b.facts['binary_length'],96);self.assertEqual(b.facts['cardinality'],4096);self.assertEqual(b.facts['minimum_distance'],32);self.assertEqual(verify_artifact(b,self.ctx).status,'exact')
 def test_searches_general_galois_ring_evaluations(self):
  ctx=RunContext(Path('.'),{'max_trials':4000,'max_rank_work':60000000,'solver_seconds':40});a=construct([],{'m':5},ctx);b=search_even_extensions([a],{'count':16,'trials':3000,'seed':0,'seconds':30,'family':'all'},ctx)
  self.assertEqual(b.facts['binary_length'],96);self.assertEqual(b.facts['cardinality'],4096);self.assertEqual(b.facts['minimum_distance'],38);self.assertEqual(verify_artifact(b,ctx).status,'exact')
 def test_discovers_translated_nonlinear_union_with_roots(self):
  a=construct([],{'m':3},self.ctx);b=find_translates([a],{'trials':20,'seed':1514,'seconds':5,'add_d_roots':True},self.ctx);self.assertEqual(b.facts['translate_layers'],8);self.assertEqual(b.facts['cardinality'],2528);self.assertEqual(verify_artifact(b,self.ctx).status,'exact')
  b.data['translates'][1]=b.data['translates'][0]
  with self.assertRaisesRegex(ValueError,'invalid binary translates'):verify_artifact(b,self.ctx)
 def test_registry(self):
  ops,_=registry();self.assertIn('galois.kerdock',ops);self.assertIn('galois.spherical',ops);self.assertIn('galois.translate_search',ops);self.assertIn('galois.evaluation_search',ops);self.assertIn('galois.find_translates',ops)
if __name__=='__main__':unittest.main()
