"""Exact 63D balanced sign refill and projected affine-quadruple trade verifier."""
from fractions import Fraction as Q
import json,sys,gzip
from pathlib import Path
from collections import defaultdict
from itertools import combinations
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from section import setup,pos

def verify():
 assert __debug__
 W,F,H,O=setup()
 J=json.loads(gzip.decompress((HERE/'witness63.json.gz').read_bytes()))
 stock=O+J
 assert len(stock)==len({s for s,a in stock})
 G=defaultdict(set)
 supported=0
 for i,(s,A) in enumerate(J):
  assert 0<s<1<<64 and s.bit_count()==9 and A and len(A)==len(set(A))
  assert all(a>=0 and a&s==a and (a&H).bit_count()*2==(s&H).bit_count() for a in A)
  assert all((a^b).bit_count()>=3 for a,b in combinations(A,2))
  assert all((s&t).bit_count()<=6 for t in W)
  for q in combinations(pos(s),4):
   t=sum(1<<j for j in q)
   if t in F:G[t].update(a&t for a in A if a&t in F[t])
  for t,T in O+J[:i]:
   p=s&t;n=p.bit_count()
   if n<=4:continue
   X={a&p for a in A};Y={b&p for b in T}
   assert not (X&Y)
   if n>=7:
    assert all(not (a^(1<<j)) in Y for a in X for j in pos(p))
   supported+=1
 a=Q(8,3)
 assert 9*a*a==64 and (9-2*3)*a*a<=32
 assert 4*a*a<=32 and 5*a*a>32
 assert 6*2*a==32 and 9*a<32 and 8*a<32
 assert 3*4*a==32 and 4*4*a>32 and (4-2)*4*a<32
 gain=sum(len(v)for _,v in J)-sum(map(len,G.values()))
 assert gain>0
 return dict(dimension=63,old=184417784,total=184417784+gain,gain=gain,
             new_supports=len(J),added_signs=sum(len(v)for _,v in J),
             removed_quadruple_signs=sum(map(len,G.values())),
             hit_old_quadruple_supports=sum(bool(v)for v in G.values()),relevant_new_support_pairs=supported)
if __name__=='__main__':print(json.dumps(verify(),indent=2))
