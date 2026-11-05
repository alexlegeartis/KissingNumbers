from pathlib import Path
from itertools import product
import unittest
from copy import deepcopy
from discovery.certification import verify_artifact
from discovery.core import RunContext
from discovery.engines.binary_translate_search import search_binary_translates, binary_translate_search
from discovery.runner import registry

class BinaryTranslateSearchTests(unittest.TestCase):
 def params(self):
  return {'generator':[[1,1,1,1]],'modulus':2,'block_map':[0,1], 'bit_width':1,
          'alphabet':[[0,1]]*4,'references':[[0,0,0,0]],'minimum_distance':2,
          'rounds':20,'oracle_steps':2,'exact_pair_limit':4,'seed':4,'seconds':5}
 def test_registry_and_exact_search(self):
  ops,kinds=registry();self.assertIn('search.binary_translates',ops);self.assertEqual(kinds['search.binary_translates'],())
  p=self.params();z=search_binary_translates(p);self.assertTrue(z['candidates']);self.assertTrue(z['rejection_witnesses'])
  C=[(0,0,0,0),(1,1,1,1)]
  for c in z['candidates']:
   self.assertTrue(c['exhaustive_cross_check']);t=c['translate']
   self.assertGreaterEqual(min(sum(x^y^w for x,y,w in zip(t,a,b)) for a,b in product(C,repeat=2)),2)
  for w in z['rejection_witnesses']:
   a,b=w['message_a'][0],w['message_b'][0]
   self.assertEqual(w['distance'],sum(x^a^b for x in w['rejected_translate']))
   self.assertLess(w['distance'],2)
 def test_heuristic_never_certifies_and_artifact_is_unverified(self):
  p=self.params();p.update(exact_pair_limit=1,rounds=5,oracle_steps=1)
  a=binary_translate_search([],p,RunContext(Path('.')))
  self.assertEqual(a.verification.status,'unverified');self.assertEqual(a.data['result']['oracle'],'heuristic')
  self.assertTrue(all(not c['exhaustive_cross_check'] for c in a.data['result']['candidates']))
 def test_receipt_replays_and_rejects_tampering(self):
  a=binary_translate_search([],self.params(),RunContext(Path('.')));r=verify_artifact(a);self.assertEqual(r.status,'unverified');self.assertGreater(r.details['witnesses_checked'],0);self.assertGreater(r.details['cross_checks_replayed'],0)
  b=deepcopy(a);b.data['result']['rejection_witnesses'][0]['distance']=99
  with self.assertRaises(ValueError):verify_artifact(b)
  b=deepcopy(a);b.data['result']['candidates'][0]['translate']=[0,0,0,0]
  with self.assertRaises(ValueError):verify_artifact(b)
 def test_refuses_bad_inputs_and_budget(self):
  p=self.params();p['alphabet'][0]=[False]
  with self.assertRaises(ValueError):search_binary_translates(p)
  p=self.params();p['rounds']=201
  with self.assertRaises(ValueError):search_binary_translates(p)
 def test_affine_ring_offsets_are_not_binary_translates(self):
  p={'generator':[[1,1,0,0],[0,1,1,0]],'modulus':3,'block_map':[0,1,3],'bit_width':2,
     'candidate_offset':[1,0,2,0],'reference_offsets':[[0,0,0,0],[1,0,1,0]],
     'alphabet':[[0,1,2,3]]*3+[[3]],'references':[[0,0,0,0],[1,2,0,0]],'minimum_distance':2,
     'rounds':4,'oracle_steps':2,'exact_pair_limit':162,'seconds':5,'seed':4}
  a=binary_translate_search([],p,RunContext(Path('.')));self.assertTrue(a.data['result']['candidates'])
  verify_artifact(a)
  G=p['generator'];m=p['block_map']
  for candidate in a.data['result']['candidates']:
   t=candidate['translate']
   for ref,off in zip(p['references'],p['reference_offsets']):
    for x,y in product(product(range(3),repeat=2),repeat=2):
     d=sum((t[j]^ref[j]^m[(sum(x[i]*G[i][j] for i in range(2))+p['candidate_offset'][j])%3]^m[(sum(y[i]*G[i][j] for i in range(2))+off[j])%3]).bit_count() for j in range(4))
     self.assertGreaterEqual(d,2)
 def test_nonadditive_map_exactly(self):
  p={'generator':[[1,1,0,0],[0,1,1,0]], 'modulus':3, 'block_map':[0,1,3], 'bit_width':2,
     'alphabet':[[0,1,2,3]]*4,'references':[[0,0,0,0]],'minimum_distance':2,
     'rounds':10,'oracle_steps':2,'exact_pair_limit':81,'seconds':5,'seed':8}
  z=search_binary_translates(p);self.assertTrue(z['candidates'])
  G=p['generator'];m=p['block_map'];C=[tuple(m[sum(a[i]*G[i][j] for i in range(2))%3] for j in range(4)) for a in product(range(3),repeat=2)]
  for c in z['candidates']:
   self.assertGreaterEqual(min(sum((t^x^y).bit_count() for t,x,y in zip(c['translate'],a,b)) for a,b in product(C,repeat=2)),2)
if __name__=='__main__':unittest.main()
