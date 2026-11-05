import math,random,unittest
from decimal import Decimal,localcontext
from pathlib import Path
from discovery.core import Artifact,RunContext
from discovery.runner import registry
from discovery.certification import verify_artifact
from discovery.engines.quadratic_group_assemblies import _positive,verify_exact,search_quadratic,certify_quadratic
class QuadraticGroupTests(unittest.TestCase):
 def setUp(self):self.ctx=RunContext(Path('.'))
 def test_exact_sign_comparison(self):
  rng=random.Random(33)
  with localcontext() as c:
   c.prec=80;root=Decimal(3).sqrt()
   for _ in range(1000):
    p=rng.randrange(-100000,100001);q=rng.randrange(-100000,100001);self.assertEqual(_positive(p,q),Decimal(p)+Decimal(q)*root>0)
 def test_exact_quadratic_boundary_and_rejection(self):
  d={'backend':'exact_quadratic_group_orbits','groups':[[],[]],'seed_a':[[1,0],[1,0]],'seed_b':[[0,0],[0,1]],'metric':[1,1]};a=Artifact('implicit_spherical_code',d);r=verify_exact(a,self.ctx);self.assertEqual(r.status,'exact');self.assertEqual(r.details['facts']['cardinality'],2);self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
  d['seed_b'][1]=[0,0]
  with self.assertRaisesRegex(ValueError,'separation'):verify_exact(Artifact('implicit_spherical_code',d),self.ctx)
 def test_random_field_pair_checks_and_large_integer_fallback(self):
  rng=random.Random(72)
  with localcontext() as c:
   c.prec=80;root=Decimal(3).sqrt()
   for _ in range(150):
    a=[[rng.randint(-3,3) for _ in range(3)] for _ in range(2)];b=[[rng.randint(-3,3) for _ in range(3)] for _ in range(2)]
    if any(not any(x+y for x,y in zip(aa,bb)) and not any(x-y for x,y in zip(aa,bb)) for aa,bb in zip(a,b)):continue
    v=[[Decimal(x)+root*Decimal(y) for x,y in zip(aa,bb)] for aa,bb in zip(a,b)];norm=[sum(x*x for x in q) for q in v]
    if min(norm)==0:continue
    dot=sum(x*y for x,y in zip(*v));bad=dot>0 and 4*dot*dot>norm[0]*norm[1]+Decimal('1e-60')
    for k in (1,10000):
     data={'backend':'exact_quadratic_group_orbits','groups':[[],[]],'seed_a':[[k*x for x in aa] for aa in a],'seed_b':[[k*x for x in bb] for bb in b],'metric':[1,1,1]}
     try:verify_exact(Artifact('implicit_spherical_code',data),self.ctx);actual=False
     except ValueError as e:
      if 'separation' not in str(e):raise
      actual=True
     self.assertEqual(actual,bad)
 def test_search_and_exact_promotion(self):
  self.assertIn('search.quadratic_group_assemblies',registry()[0]);p={'dimension':2,'orbit_specs':[{'generators':[]},{'generators':[]}],'initial_a':[[1,0],[1,0]],'initial_b':[[0,0],[0,1]],'keep_initial':True,'attempts':1,'max_iterations':1,'seconds':4,'scale':100}
  r=search_quadratic([],p,self.ctx);self.assertEqual(r.verification.status,'unverified');a=certify_quadratic([r],{},self.ctx);self.assertEqual(a.facts['cardinality'],2)
if __name__=='__main__':unittest.main()
