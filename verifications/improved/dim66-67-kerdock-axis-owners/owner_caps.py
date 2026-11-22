"""Exact Delsarte certificates for the fixed doubly-even dense owner scopes."""
from fractions import Fraction as Q
from math import comb
C = lambda n, k: comb(n, k) if 0 <= k <= n else 0
def K(j,w):
    return Q(sum((-1)**i*C(w,i)*C(64-w,j-i) for i in range(j+1)),C(64,j))
cases = [
 (24,131072,{1:Q(1413,40),2:Q(591,2),3:Q(333901,120),4:Q(9455),
              5:Q(4290679,120),6:Q(111569,2),7:Q(3235501,120)}),
 (28,4096,{1:Q(4234169,82322),2:Q(156),3:Q(853187549,535093),
           4:Q(1891),9:Q(430321633,1070186)})]
for d,N,y in cases:
    F=lambda w: 1+sum(x*K(j,w)for j,x in y.items())
    assert all(x>0 for x in y.values())
    assert F(0)==N and all(F(w)<=0 for w in range(d,65,4))
    print(f"all-multiples-of-four binary distance set, minimum {d}: upper bound {N}")
