from itertools import combinations
from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.spectral_graph_search import search
from discovery.runner import registry

class SpectralGraphSearchTests(unittest.TestCase):
 def setUp(self):self.ctx=RunContext(Path('.'))
 def test_triangular_graph_representation(self):
  edges=list(combinations(range(8),2));A=[[int(i!=j and len(set(edges[i])&set(edges[j]))==1) for j in range(28)] for i in range(28)]
  r=search([] ,{'adjacency':A,'max_dimension':7,'max_sectors':1,'include_antipodes':True,'seconds':5},self.ctx)
  self.assertEqual(r.facts['cardinality'],56);self.assertEqual(r.facts['rank_exact'],7);self.assertEqual(r.facts['max_inner_product'],[1,3]);self.assertEqual(verify_artifact(r,self.ctx).status,'exact')
 def test_simplex_and_antipodes(self):
  A=[[int(i!=j) for j in range(4)] for i in range(4)]
  r=search([] ,{'adjacency':A,'max_dimension':3,'max_sectors':1,'seconds':5},self.ctx);self.assertEqual(r.facts['cardinality'],8);self.assertEqual(r.facts['rank_exact'],3)
 def test_quadratic_pentagon_projector(self):
  A=[[int((i-j)%5 in (1,4)) for j in range(5)] for i in range(5)]
  r=search([],{'adjacency':A,'max_dimension':2,'max_sectors':1,'seconds':5},self.ctx)
  self.assertEqual(r.kind,'quadratic_projector_code');self.assertEqual(r.facts['cardinality'],5);self.assertEqual(r.facts['ambient_dimension'],2);self.assertEqual(verify_artifact(r,self.ctx).status,'exact')
 def test_icosahedron_quadratic_graph_recovers_12(self):
  A=[[0]*12 for _ in range(12)]
  def edge(i,j):A[i][j]=A[j][i]=1
  for i in range(5):
   edge(0,1+i);edge(11,6+i);edge(1+i,1+(i+1)%5);edge(6+i,6+(i+1)%5);edge(1+i,6+i);edge(1+i,6+(i-1)%5)
  r=search([],{'adjacency':A,'max_dimension':3,'max_sectors':1,'seconds':10},RunContext(Path('.'),{'solver_seconds':15}))
  self.assertEqual(r.kind,'quadratic_projector_code');self.assertEqual(r.facts['cardinality'],12);self.assertEqual(r.facts['ambient_dimension'],3);self.assertEqual(verify_artifact(r,self.ctx).status,'exact')
 def test_rejects_unsupported_spectrum_and_asymmetry(self):
  A=[[int((i-j)%7 in (1,6)) for j in range(7)] for i in range(7)]
  with self.assertRaisesRegex(ValueError,'quadratic'):search([],{'adjacency':A,'seconds':2},self.ctx)
  with self.assertRaisesRegex(ValueError,'symmetric'):search([],{'adjacency':[[0,1],[0,0]]},self.ctx)
 def test_rejects_nonprojector(self):
  from discovery.core import Artifact
  from discovery.engines.quadratic_projector import verify_artifact as check
  a=Artifact('quadratic_projector_code',{'radicand':5,'projection':[[[[1,2],0]]],'selected':[[0,1]]})
  with self.assertRaisesRegex(ValueError,'not an exact projector'):check(a,self.ctx)
 def test_registry(self):self.assertIn('search.spectral_graph',registry()[0])
if __name__=='__main__':unittest.main()
