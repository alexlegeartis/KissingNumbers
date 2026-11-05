from copy import deepcopy
from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.ring_orbit_search import search_ring_orbits,ring_orbit_search
from discovery.runner import registry

class RingOrbitSearchTests(unittest.TestCase):
 def test_perfect_matching_is_independently_replayed(self):
  p={'generator':[[1]],'alpha':[1],'shift':[0],'minimum_distance':1,'max_orbits':2,'attempts':1,'steps':1,'seed':1,'seconds':4}
  ops,kinds=registry();self.assertIn('search.ring_orbit_assemblies',ops);self.assertEqual(kinds['search.ring_orbit_assemblies'],())
  a=ring_orbit_search([],p,RunContext(Path('.')));r=a.data['result'];self.assertEqual(len(r['matching']),2);self.assertEqual(r['whole_orbit_count_upper_bound'],2);self.assertEqual(r['orbit_size'],4)
  v=verify_artifact(a);self.assertEqual(v.status,'unverified');self.assertEqual(v.details['matching_edges_checked'],2)
  b=deepcopy(a);b.data['result']['witnesses'][0]['distance']=1
  with self.assertRaises(ValueError):verify_artifact(b)
  b=deepcopy(a);b.data['result']['matching'].append(b.data['result']['matching'][0])
  with self.assertRaises(ValueError):verify_artifact(b)
 def test_nontrivial_octonary_scaling_symmetry(self):
  p={'generator':[[1,1],[0,1]],'alpha':[5,1],'shift':[0,0],'minimum_distance':1,'max_orbits':4,'attempts':1,'steps':1,'seed':3,'seconds':4}
  a=ring_orbit_search([],p,RunContext(Path('.')));r=a.data['result'];self.assertEqual(r['subgroup_symmetries'],4);self.assertEqual(len(r['matching']),4);self.assertEqual(verify_artifact(a).details['matching_edges_checked'],4)
 def test_provisional_set_is_independent_in_observed_graph(self):
  from discovery.engines.ring_orbit_search import _parameters,_expand
  p={'generator':[[1,0,1],[0,1,1]],'alpha':[5,1,5],'shift':[1,0,2],'minimum_distance':3,'max_orbits':4,'attempts':2,'steps':20,'seed':8,'seconds':4}
  ctx=RunContext(Path('.'));r=search_ring_orbits(p,ctx);G,alpha,shift,target,bits,K=_parameters(p,ctx);edges=_expand(r['witnesses'],G,alpha,shift,target,bits,K,ctx)
  A=set(r['provisional_candidate_orbits']);B=set(r['provisional_source_orbits']);self.assertFalse(any(a in A and b in B for a,b in edges));self.assertEqual(len(A)+len(B),r['whole_orbit_count_upper_bound'])
 def test_unobserved_conflicts_are_not_certified(self):
  p={'generator':[[1,1]],'alpha':[1,1],'shift':[0,1],'minimum_distance':1,'max_orbits':2,'attempts':1,'steps':1,'seed':7,'seconds':4}
  r=search_ring_orbits(p);self.assertEqual(r['matching'],[]);self.assertEqual(r['whole_orbit_count_upper_bound'],4)
  self.assertIn('provisional',r['scope'].lower());self.assertEqual(len(r['provisional_candidate_orbits'])+len(r['provisional_source_orbits']),4)
 def test_seed_witnesses_are_checked_and_used(self):
  p={'generator':[[1]],'alpha':[1],'shift':[0],'minimum_distance':1,'max_orbits':1,'attempts':1,'steps':1,'seed':1,'seconds':4,'initial_witnesses':[{'candidate':0,'source_coefficients':[0],'distance':0}]}
  a=ring_orbit_search([],p,RunContext(Path('.')));self.assertEqual(len(a.data['result']['matching']),2);self.assertEqual(verify_artifact(a).status,'unverified')
  p['initial_witnesses'][0]['distance']=1
  with self.assertRaises(ValueError):search_ring_orbits(p)
 def test_freeness_and_units_validated(self):
  p={'generator':[[1,1],[1,1]],'alpha':[1,1],'shift':[0,0],'minimum_distance':2}
  with self.assertRaises(ValueError):search_ring_orbits(p)
  p={'generator':[[1]],'alpha':[3],'shift':[0],'minimum_distance':2}
  with self.assertRaises(ValueError):search_ring_orbits(p)
if __name__=='__main__':unittest.main()
