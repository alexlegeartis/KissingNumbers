from pathlib import Path
import unittest
from itertools import combinations
from discovery.engines.signed_exchange import _bad
from discovery.core import RunContext
from discovery.certification import verify_artifact
from discovery.engines.layers import signed_supports
from discovery.engines.signed_exchange import signed_exchange, _lines
from discovery.runner import registry

class SignedExchangeTests(unittest.TestCase):
 def setUp(self):
  self.ctx=RunContext(Path('.'))
  self.cert={'format':'signed_support_layers_v1','dimension':32,'dense_length':32,
   'dense_generator':[(1<<32)-1],'full_even_supports':[255],
   'exceptional_antipodal_lines':[],'expected_cardinality':2114}
  self.base=signed_supports([] ,{'certificate':self.cert},self.ctx)
 def test_constructs_and_verifies_new_line_layer(self):
  ops,kinds=registry();self.assertIn('search.signed_fiber_exchanges',ops);self.assertEqual(kinds['search.signed_fiber_exchanges'],('implicit_spherical_code',))
  candidate=[list(x) for x in _lines(65280,0)]
  out=signed_exchange([self.base],{'candidate_lines':candidate,'seconds':4},self.ctx)
  self.assertEqual(out.facts['cardinality'],2242);self.assertEqual(verify_artifact(out,self.ctx).status,'exact')
 def test_neutral_odd_fiber_does_not_claim_gain(self):
  odd=[list(x) for x in _lines(255,1)]
  out=signed_exchange([self.base],{'candidate_lines':odd,'flexible_supports':[255],'seconds':4},self.ctx)
  self.assertEqual(out.facts['cardinality'],2114);self.assertEqual(verify_artifact(out,self.ctx).status,'exact');self.assertEqual(out.data['provenance']['classification'],'signed_fiber_search');self.assertTrue(out.data['provenance']['solver_status'])
 def test_joint_parity_model_matches_brute_force(self):
  candidate=_lines(0xff00,0)[:4]
  candidate+=_lines(0x033f,1)[:4]
  full=[_lines(255,0),_lines(255,1)];best=0
  for mask in range(1<<len(candidate)):
   selected=[candidate[i] for i in range(len(candidate)) if mask>>i&1]
   if any(_bad(a,b) for a,b in combinations(selected,2)):continue
   cost=min(sum(any(_bad(x,y) for x in selected) for y in fiber) for fiber in full)
   best=max(best,len(selected)-cost)
  out=signed_exchange([self.base],{'candidate_lines':[list(x) for x in candidate], 'flexible_supports':[255],'seconds':4},self.ctx)
  self.assertEqual(out.facts['cardinality'],2114+2*best);self.assertEqual(verify_artifact(out,self.ctx).status,'exact')
 def test_can_return_neutral_switched_fiber(self):
  odd=[list(x) for x in _lines(255,1)]
  out=signed_exchange([self.base],{'candidate_lines':odd,'flexible_supports':[255],'min_new_lines':64,'return_neutral':True,'seconds':4},self.ctx)
  self.assertEqual(out.facts['cardinality'],2114);self.assertEqual(out.data['certificate']['full_even_supports'],[])
  self.assertEqual(len(out.data['certificate']['exceptional_antipodal_lines']),64);self.assertEqual(verify_artifact(out,self.ctx).status,'exact')
 def test_explicit_triangle_covers_for_all_overlap_types(self):
  old=255
  for r in (5,6,7):
   common=list(range(r));old_extra=list(range(r,8));new_extra=list(range(8,8+8-r));new=sum(1<<i for i in common+new_extra)
   for old_parity in (0,1):
    oldlines=_lines(old,old_parity);seen=[set(),set()]
    for _,v in oldlines:
     candidate=[]
     for parity in (0,1):
      w=sum(1<<j for j in common if v>>j&1)
      for a,b in zip(old_extra[:-1],new_extra[:-1]):
       if v>>a&1:w|=1<<b
      if (w.bit_count()&1)!=parity:w|=1<<new_extra[-1]
      candidate.append((new,w));seen[parity].add(w)
     self.assertTrue(_bad((old,v),candidate[0]));self.assertTrue(_bad((old,v),candidate[1]));self.assertTrue(_bad(candidate[0],candidate[1]))
    self.assertEqual([len(x) for x in seen],[64,64])
 def test_explicit_quadruple_covers_for_overlapping_new_supports(self):
  first=255
  for r in (5,6,7):
   common=list(range(r));extra1=list(range(r,8));extra2=list(range(8,8+8-r));second=sum(1<<i for i in common+extra2);seen=[set() for _ in range(4)]
   for _,v in _lines(first,0):
    lines=[]
    for support,extra in ((first,extra1),(second,extra2)):
     for parity in (0,1):
      w=sum(1<<j for j in common if v>>j&1)
      for a,b in zip(extra1[:-1],extra[:-1]):
       if v>>a&1:w|=1<<b
      if (w.bit_count()&1)!=parity:w|=1<<extra[-1]
      seen[len(lines)].add(w);lines.append((support,w))
    self.assertTrue(all(_bad(a,b) for a,b in combinations(lines,2)))
   self.assertEqual([len(x) for x in seen],[64]*4)
 def test_two_support_cuts_match_exhaustive_search(self):
  import random
  rng=random.Random(314159);pool=_lines(0x033f)+_lines(0x053f);old=[_lines(255,0),_lines(255,1)]
  for trial in range(5):
   candidate=rng.sample(pool,8);best=0
   for mask in range(1<<len(candidate)):
    selected=[candidate[i] for i in range(len(candidate)) if mask>>i&1]
    if any(_bad(a,b) for a,b in combinations(selected,2)):continue
    cost=min(sum(any(_bad(x,y) for x in selected) for y in fiber) for fiber in old)
    best=max(best,len(selected)-cost)
   out=signed_exchange([self.base],{'candidate_lines':[list(x) for x in candidate], 'flexible_supports':[255],'seconds':4},self.ctx)
   self.assertEqual(out.facts['cardinality'],2114+2*best)
 def test_nine_coordinate_cut_matches_exhaustive_search(self):
  import random
  rng=random.Random(919);old=[_lines(255,0),_lines(255,1)];pool=[]
  for missing in (0,1,2):pool+=_lines(((1<<9)-1)^(1<<missing))
  for trial in range(5):
   candidate=rng.sample(pool,8);best=0
   for mask in range(1<<len(candidate)):
    selected=[candidate[i] for i in range(len(candidate)) if mask>>i&1]
    if any(_bad(a,b) for a,b in combinations(selected,2)):continue
    cost=min(sum(any(_bad(x,y) for x in selected) for y in fiber) for fiber in old)
    best=max(best,len(selected)-cost)
   out=signed_exchange([self.base],{'candidate_lines':[list(x) for x in candidate], 'flexible_supports':[255],'seconds':4},self.ctx)
   self.assertEqual(out.facts['cardinality'],2114+2*best)
 def test_common_six_bound_matches_exhaustive_search(self):
  import random
  rng=random.Random(669);old=[_lines(255,0),_lines(255,1)];pool=[]
  for extra in ((6,8),(7,9),(8,9)):
   support=sum(1<<j for j in list(range(6))+list(extra));pool+=_lines(support)
  for trial in range(5):
   candidate=rng.sample(pool,8);best=0
   for mask in range(1<<len(candidate)):
    selected=[candidate[i] for i in range(len(candidate)) if mask>>i&1]
    if any(_bad(a,b) for a,b in combinations(selected,2)):continue
    cost=min(sum(any(_bad(x,y) for x in selected) for y in fiber) for fiber in old)
    best=max(best,len(selected)-cost)
   out=signed_exchange([self.base],{'candidate_lines':[list(x) for x in candidate], 'flexible_supports':[255],'seconds':4},self.ctx)
   self.assertEqual(out.facts['cardinality'],2114+2*best)
 def test_matching_bound_prunes_finite_pool_exactly(self):
  c=dict(self.cert);c['full_even_supports']=[255,65280];c['expected_cardinality']=2242
  a=signed_supports([] ,{'certificate':c},self.ctx)
  s1=0x033f;s2=(0x3f<<8)|(3<<16);candidate=[list(_lines(s1,0)[0]),list(_lines(s2,0)[0])]
  out=signed_exchange([a],{'candidate_lines':candidate,'flexible_supports':[255,65280],'seconds':4},self.ctx)
  self.assertEqual(out.facts['cardinality'],2242);self.assertIn('structural bound',out.data['provenance']['solver_status']);self.assertEqual(verify_artifact(out,self.ctx).status,'exact')
 def test_invalid_flexible_and_candidate_signs_rejected(self):
  with self.assertRaises(ValueError):signed_exchange([self.base],{'candidate_supports':[65280],'flexible_supports':[17]},self.ctx)
  with self.assertRaises(ValueError):signed_exchange([self.base],{'candidate_lines':[[65280,1]]},self.ctx)
if __name__=='__main__':unittest.main()
