"""Exact P48 11-design certificates for particular 47D layer families."""
from fractions import Fraction as F
from pathlib import Path
import itertools,json,math
import sympy as S
N=52416000;dim=48

def df(n):return math.prod(range(n,0,-2)) if n>0 else 1
def moment(a,b,g11,g22,g12):
 if (a+b)%2:return F(0)
 k=(a+b)//2
 z=sum(math.comb(a,t)*math.comb(b,t)*math.factorial(t)*df(a-t-1)*df(b-t-1)*g11**((a-t)//2)*g22**((b-t)//2)*g12**t for t in range(min(a,b)+1) if not (a-t)%2 and not (b-t)%2)
 return F(6**k*z,math.prod(dim+2*i for i in range(k)))
# u and v are distinct minimal vectors of norm six with u.v=-2.
# Four exceptional vectors +-u,+-v are removed, leaving s,t in [-3,3].
c=-2;pos=list(itertools.product(range(-3,4),repeat=2));exps=[(a,b) for a in range(7) for b in range(7) if a+b<=11]
A=S.Matrix([[s**a*t**b for s,t in pos] for a,b in exps]);B=S.Matrix([S.Rational(N)*S.Rational(moment(a,b,6,6,c))-(6**a*c**b+(-6)**a*(-c)**b+c**a*6**b+(-c)**a*(-6)**b) for a,b in exps]);v,p=A.gauss_jordan_solve(B);assert p.rows==1;tau=p[0];assert v[pos.index((3,3))]==tau and v[pos.index((3,2))]==-6*tau
v=[int(x.subs(tau,0)) for x in v];assert min(v)>=0
counts={r:sum(v[i] for i,(s,t) in enumerate(pos) if s+t==r) for r in range(-6,7)}
for s,t in [(6,c),(-6,-c),(c,6),(-c,-6)]:counts[s+t]=counts.get(s+t,0)+1
assert sum(counts.values())==N and counts[0]==20582240 and counts[-1]==counts[1]==12816384
assert all(counts[r]==0 for r in (-6,-5,5,6))
print('normal norm 8: central section',counts[0],'each projected +/-1 layer',counts[1])
# For a minimal x in layer c0, classify other shell vectors by s=x.y and
# t=(u+v).y. Outside +-x, s in [-3,3], and all t lie in [-4,4].
pos2=list(itertools.product(range(-3,4),range(-4,5)));duals=json.loads((Path(__file__).parent/'duals.json').read_text());proved={}
for c0,target in [(0,1),(1,0),(1,-1)]:
 for sign in (1,-1):
  coeff=[(int(a),int(b),F(y)) for a,b,y in duals[f'{c0},{target},{sign}']]
  assert all(a+b<=11 for a,b,_ in coeff)
  for s,t in pos2:
   lhs=sum(y*F(s**a,3**a)*F(t**b,4**b) for a,b,y in coeff)
   assert lhs <= sign*int((s,t)==(3,target)),(c0,target,sign,s,t)
  total=F(0)
  for a,b,y in coeff:
   excluded=F(6**a*c0**b+(-6)**a*(-c0)**b,N)
   total += y*(moment(a,b,6,8,c0)-excluded)/3**a/4**b
  proved[c0,target,sign]=total*N
expected={(0,1,1):F(8448),(0,1,-1):F(-12884),(1,0,1):F(13113),(1,0,-1):F(-29751,2),(1,-1,1):F(3666),(1,-1,-1):F(-14053,3)}
assert proved==expected
print('exact cross degrees: central to cap 8448..12884; cap to central 13113..14875; opposite caps 3666..4684')
# The edge condition for the +/-1 projected layers is original product 3.
# By Konig, a min vertex cover has at least |E|/Delta vertices.
n=counts[1];matching=(3666*n+4683)//4684;bound=2*n-matching
assert matching==10030928 and bound==15601840 and bound<23766960
print('both projected +/-1 layers at most',bound,'points')
# Degree inequalities for any independent subset of 0,+1,-1 layers.
# Exact LP dual: first central-cap inequality twice and cap-cap once.
# With equal cap weights, cap occupancy <= 4684*n/(4684+3666).
# a0 <= (14875/8448)*(n-a+), and symmetric for a-.
pmax=F(4684*n,4684+3666);upper=F(14875,8448)*n+(2-F(14875,8448))*pmax
assert upper>23766960 and math.floor(upper)==24286647
print('0,+/-1 degree relaxation upper',math.floor(upper),'(not a construction)')
# Separately, use a minimal normal u. Its orthogonal section is the known
# 23,766,960-point 47D configuration. Project the u.x=+1 layer.
# Any added projected point is blocked exactly by source z with x.z=3.
for c in (0,1):
 B=S.Matrix([S.Rational(N)*S.Rational(moment(a,b,6,6,c))-(6**a*c**b+(-6)**a*(-c)**b+c**a*6**b+(-c)**a*(-6)**b) for a,b in exps]);v,p=A.gauss_jordan_solve(B);assert p.rows==1;tau=p[0];lo=S.Rational(0);hi=S.oo
 for f in v:
  slope=S.diff(f,tau);offset=f.subs(tau,0)
  if slope>0:lo=max(lo,-offset/slope)
  elif slope<0:hi=min(hi,offset/(-slope))
 lo,hi=int(S.ceiling(lo)),int(S.floor(hi))
 if c==0:
  assert (lo,hi)==(0,96) and v[pos.index((1,3))]==8280+15*tau
 else:
  assert (lo,hi)==(49,57) and v[pos.index((0,3))]==16660-20*tau
# Each new cap point has >=15,520 old blockers; each old point blocks
# <=9,720 such cap points. Shared-deletion counting prohibits a net gain.
assert 16660-20*57==15520 and 8280+15*96==9720 and 15520>9720
print('minimal-normal section: a single +/-1 projected cap alone cannot improve the 23,766,960 source')
