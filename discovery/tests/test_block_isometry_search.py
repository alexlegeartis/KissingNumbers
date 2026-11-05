from copy import deepcopy
from itertools import product
from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.block_isometry_search import search_block_isometries, block_isometry_search
from discovery.runner import registry

class BlockIsometryTests(unittest.TestCase):
 def params(self):
  return {'generator':[[1,1,1,1]],'modulus':2,'bit_width':1,'block_map':[0,1],
          'map_menu':[[0,1],[1,0]],'groups':[0,1,2,3],'references':[[0,0,0,0]],
          'minimum_distance':2,'rounds':20,'oracle_steps':2,'exact_pair_limit':4,'seed':4,'seconds':5}
 def test_exact_search_and_receipt(self):
  ops,kinds=registry();self.assertIn('search.block_isometries',ops);self.assertEqual(kinds['search.block_isometries'],())
  a=block_isometry_search([],self.params(),RunContext(Path('.')));r=a.data['result'];self.assertTrue(r['candidates']);self.assertTrue(r['rejection_witnesses'])
  check=verify_artifact(a);self.assertEqual(check.status,'unverified');self.assertGreater(check.details['cross_checks_replayed'],0)
  for c in r['candidates']:
   choice=c['map_choices'];self.assertGreaterEqual(min(sum(choice[j]^x^y for j in range(4)) for x,y in product(range(2),repeat=2)),2)
  b=deepcopy(a);b.data['result']['rejection_witnesses'][0]['distance']=9
  with self.assertRaises(ValueError):verify_artifact(b)
 def test_checks_internal_cell_isometry(self):
  p=self.params();p['map_menu']=[[0,0]]
  with self.assertRaises(ValueError):search_block_isometries(p)
  p=self.params();p['groups']=[0,1,2,4]
  with self.assertRaises(ValueError):search_block_isometries(p)
 def test_affine_offsets_and_group_ties(self):
  p={'generator':[[1,1,0]],'modulus':3,'bit_width':2,'block_map':[0,1,3],
     'candidate_offset':[1,0,0],'reference_offsets':[[2,1,0],[0,1,0]],
     'map_menu':[[0,1,3],[3,2,0]],'groups':[0,0,1],
     'references':[[0,0,0],[1,1,0]],'minimum_distance':1,'rounds':4,
     'oracle_steps':2,'exact_pair_limit':18,'seed':8,'seconds':5}
  a=block_isometry_search([],p,RunContext(Path('.')));self.assertTrue(a.data['result']['candidates']);verify_artifact(a)
  for c in a.data['result']['candidates']:
   choices=c['map_choices'];
   for ri,ref in enumerate(p['references']):
    for x,y in product(range(3),repeat=2):
     d=sum((p['map_menu'][choices[p['groups'][j]]][(x*p['generator'][0][j]+p['candidate_offset'][j])%3]^ref[j]^p['block_map'][(y*p['generator'][0][j]+p['reference_offsets'][ri][j])%3]).bit_count() for j in range(3))
     self.assertGreaterEqual(d,1)
 def test_exact_information_set_overlap_rejection(self):
  p=self.params();p.update(map_menu=[[0,1]],groups=[0,0,0,0],exact_pair_limit=1,rounds=1,oracle_steps=1)
  r=search_block_isometries(p);self.assertEqual(len(r['rejection_witnesses']),1);self.assertEqual(r['rejection_witnesses'][0]['distance'],0);self.assertFalse(r['candidates'])
 def test_anchored_search_finds_noncident_blocker(self):
  p=self.params();p.update(map_menu=[[0,1]],groups=[0,0,0,0],references=[[1,0,0,0]],exact_pair_limit=1,rounds=1,oracle_steps=100)
  r=search_block_isometries(p);self.assertEqual(len(r['rejection_witnesses']),1);w=r['rejection_witnesses'][0];self.assertEqual(w['distance'],1);self.assertEqual(w['message_a'],[0])
 def test_heuristic_does_not_certify(self):
  p=self.params();p.update(exact_pair_limit=1,oracle_steps=1,rounds=3)
  z=search_block_isometries(p);self.assertEqual(z['oracle'],'heuristic');self.assertTrue(all(not x['exhaustive_cross_check'] for x in z['candidates']))
if __name__=='__main__':unittest.main()
