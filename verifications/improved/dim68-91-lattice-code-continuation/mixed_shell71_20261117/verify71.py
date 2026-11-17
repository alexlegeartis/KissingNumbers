from fractions import Fraction as F
from itertools import combinations_with_replacement
import json, math
N8=6218175600;N10=15281788354560
m8={0:2603658750,1:1512243200,2:280928256,3:13959168,4:127800,8:1}
m10={0:5723204788224,1:3708252979200,2:975222374400,3:93206937600,4:2595532800,5:13959168}
for s,N,M in [(8,N8,m8),(10,N10,m10)]:
 for k in range(6):
  observed=M[0] if k==0 else 0
  observed+=2*sum(c*j**(2*k) for j,c in M.items() if j)
  expected=F(N*(8*s)**k*math.prod(range(1,2*k,2)),math.prod(range(72,72+2*k,2)))
  assert observed==expected,(s,k,observed,expected)
slices=[(8,3,m8[3]),(8,4,m8[4]),(10,4,m10[4]),(10,5,m10[5])]
checks=[]
for (s,j,_),(t,k,_) in combinations_with_replacement(slices,2):
 ri=F(s)-F(j*j,8);rj=F(t)-F(k*k,8)
 bound=F(s+t-8,2)-F(j*k,8)
 assert ri>0 and rj>0
 assert bound<=0 or 4*bound*bound<=ri*rj
 checks.append({'types':[[s,j],[t,k]],'residual_norms':[str(ri),str(rj)],'dot_upper':str(bound),'squared_margin':str(ri*rj-4*bound*bound)})
total=sum(c for _,_,c in slices)
assert total==2623578936
print(json.dumps({'status':'PASS','dimension':71,'cardinality':total,'gain_over_repository':total-2603658750,'checks':checks,'premises':['Existence of extremal even unimodular rank72 lattice with minimum8','Both norm8 and norm10 shells are spherical11designs','Theta-series coefficients N8,N10 as stated']},indent=2))
