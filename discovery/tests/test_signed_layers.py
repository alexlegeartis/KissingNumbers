from collections import Counter
from copy import deepcopy
from decimal import Decimal,localcontext
from itertools import product
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
from discovery.core import Artifact,RunContext,Verification
from discovery.engines import signed_layers as layers


def qcode(A,B=None):
 return Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':A,'radical_coordinates':B or [[0]*len(A[0]) for _ in A]})


def classical6():
 rows=[]
 for supp in ((0,1,2,3),(0,1,4,5),(2,3,4,5)):
  for signs in product((-1,1),repeat=4):
   row=[0]*6
   for i,s in zip(supp,signs):row[i]=s
   rows.append(row)
 return Artifact('spherical_code',{'coordinates':rows,'norm_squared':4})


class SignedLayerTests(unittest.TestCase):
 def setUp(self):self.context=RunContext(Path('.'))
 def test_exact_quadratic_sign_against_high_precision(self):
  with localcontext() as ctx:
   ctx.prec=60;root=Decimal(3).sqrt()
   for a in range(-15,16):
    for b in range(-15,16):
     value=Decimal(a)+Decimal(b)*root
     self.assertEqual(layers.sign_sqrt3(a,b),(value>0)-(value<0))
 def test_variable_norms_and_nonantipodal_points(self):
  a=qcode([[1,0],[-100,0],[0,3]])
  facts=layers.verify_artifact(a).details['facts']
  self.assertEqual(facts['cardinality'],3);self.assertEqual(facts['rank_exact'],2)
  self.assertNotIn('norm_squared',facts)
 def test_sqrt3_cap_geometry_and_duplicate_rejection(self):
  a=qcode([[0,2],[0,-2],[0,0]],[[0,0],[0,0],[2,0]])
  self.assertEqual(layers.verify_artifact(a).status,'exact')
  for bad in (qcode([[1,0],[2,0]]),qcode([[1,0],[1,1]]),qcode([[0,0]]),qcode([[1,0]],[[1,0]])):
   with self.assertRaises(ValueError):layers.verify_artifact(bad)
 def test_h2_switch_is_actual_nonclassical_signed_code(self):
  result=layers.conditional_swap([classical6()],{'swap':[2,5],'conditions':[{'0':1}],'support_masks':[15,51]},self.context)
  counts=Counter(sum(1<<i for i,x in enumerate(row) if x) for row in result.data['coordinates'])
  self.assertEqual(sorted(counts.values()),[8,8,8,8,16]);self.assertEqual(result.facts['cardinality'],48)
  self.assertNotEqual(set(map(tuple,result.data['coordinates'])),set(map(tuple,classical6().data['coordinates'])))
 def test_coupled_lift_is_dimension_parametrized_and_checked(self):
  source=Artifact('spherical_code',{'coordinates':[[1],[-1]],'norm_squared':1})
  result=layers.coupled_lift([source],{'equator_scale':2,'cap_rational':[],'cap_radical':[],'heights':[-2,2],'pole_amplitude':2},self.context)
  self.assertEqual(result.facts['cardinality'],4);self.assertEqual(result.facts['rank_exact'],2)
 def test_free_zone_join_replays_classified_816_equator(self):
  second=layers.conditional_swap([classical6()],{'swap':[2,5],'conditions':[{'0':1}],'support_masks':[15,51]},self.context)
  matchings=[[[5,c],[(c+1)%5,(c-1)%5],[(c+2)%5,(c-2)%5]] for c in range(5)]
  joined=layers.free_zone_join([classical6(),second],{'left_matchings':matchings,'right_matchings':matchings},self.context)
  self.assertEqual(joined.facts['cardinality'],816)
  self.assertEqual(joined.facts['rank_exact'],12)
  bad=[[[0,1],[0,1]]]
  with self.assertRaises(ValueError):layers.free_zone_join([classical6(),second],{'left_matchings':bad,'right_matchings':bad},self.context)
 def test_coupled_output_preflight_precedes_coordinate_allocation(self):
  source=Artifact('spherical_code',{'coordinates':[[1],[-1]],'norm_squared':1})
  params={'equator_scale':2,'cap_rational':[],'cap_radical':[],'heights':[-2,2],'pole_amplitude':2}
  for budget in ({'max_dimension':1},{'max_pairs':2},{'max_rational_bits':1}):
   with self.subTest(budget=budget),patch.object(layers,'verify_artifact',side_effect=AssertionError('output verifier reached')):
    with self.assertRaises(ValueError):layers.coupled_lift([source],params,RunContext(Path('.'),budget))
 def test_real_1154_nonclassical_seed_exact(self):
  path=Path(__file__).resolve().parents[2]/'research/coupled13/seed_0_2.npz'
  with np.load(path) as data:a=qcode(data['A'].tolist(),data['B'].tolist())
  report=layers.verify_artifact(a)
  self.assertEqual((report.details['facts']['cardinality'],report.details['facts']['rank_exact']),(1154,13))
  self.assertEqual(report.details['pairs_checked'],665281)
  self.assertEqual(report.details['norm_histogram'],{'16':1154})
 def test_forged_facts_and_limits(self):
  a=qcode([[1,0],[-1,0]]);a.facts={'cardinality':999};a.verification=Verification('exact','forged')
  with self.assertRaises(ValueError):layers.verify_artifact(a)
  with patch.object(layers.np,'array',side_effect=AssertionError('allocation reached')):
   with self.assertRaises(ValueError):layers.verify_artifact(qcode([[1,0],[-1,0]]),RunContext(Path('.'),{'max_dimension':1}))
 def test_malformed_switch_and_quadratic_data(self):
  for change in ('boolean','radicand','ragged','float'):
   a=qcode([[1,0],[-1,0]])
   if change=='boolean':a.data['rational_coordinates'][0][0]=True
   elif change=='radicand':a.data['radicand']=2
   elif change=='ragged':a.data['radical_coordinates'][1]=[0]
   else:a.data['rational_coordinates'][0][0]=1.0
   with self.subTest(change=change),self.assertRaises(ValueError):layers.verify_artifact(a)
  with self.assertRaises(ValueError):layers.conditional_swap([classical6()],{'swap':[2,2],'conditions':[{'0':1}]},self.context)


if __name__=='__main__':unittest.main()
