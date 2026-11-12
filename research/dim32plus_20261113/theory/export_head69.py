"""Freeze exact rational geometry and a population dual for the 69D candidate."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import factorial, prod, gcd, lcm
import numpy as np
import sympy as sp
from scipy.optimize import linprog
import json

p=Path(__file__).parent;root=p.resolve().parents[2]
base=root/'verifications/improved/dim68-69-gamma72-cross-sections'
W=json.loads((base/'weyl_averaged/certificate4.json').read_text())
T=json.loads((base/'affine68_69/certificate4.json').read_text())
C=json.loads((p/'heads_k4_result.json').read_text())
G=sp.Matrix(W['gram']);H=G.inv();N=6218175600;k=4
offsets=[a['profile'] for a in C['slices']]
heights=list(map(F,['1/8','-1','-2/15','0','1/2','-1/2','-4/15','4/15','1/2','1']))
q=lambda a,b:F((sp.Matrix(a).T*H*sp.Matrix(b))[0])
alpha2=[(8-h*h)/(8-q(c,c)) for c,h in zip(offsets,heights)]
gaps=[]
for i,c in enumerate(offsets):
 for j,d in enumerate(offsets[:i+1]):
  R=4-heights[i]*heights[j];P=4-q(c,d)
  gap=R*R-alpha2[i]*alpha2[j]*P*P
  assert R>=0 and (P<=0 or gap>=0),(i,j,gap)
  gaps.append({'i':i,'j':j,'squared_gap':str(gap)})
print('EXACT GEOMETRY PASS',len(gaps),flush=True)
cells=list(map(tuple,W['cells']));groups=W['groups']
roots=[tuple(G*sp.Matrix(a)) for a in product(range(-2,3),repeat=4) if (sp.Matrix(a).T*G*sp.Matrix(a))[0]==8]
x=sp.symbols('x:4');quad=sum(G[i,j]*x[i]*x[j] for i in range(4) for j in range(4))
polys=[sp.Poly(quad**m,*x) for m in range(6)]
M={0:2603658750,1:1512243200,2:280928256,3:13959168,4:127800}
rows=[];rhs=[]
for tag in T['tags']:
 if tag[0]=='mom':
  powers=tag[1];degree=sum(powers);m=degree//2
  val=lambda c:prod(int(v)**a for v,a in zip(c,powers))
  coeff=int(polys[m].coeff_monomial(prod(z**a for z,a in zip(x,powers))))
  multi=factorial(degree)//prod(factorial(a) for a in powers)
  b=F(N*8**m*prod(range(1,2*m,2))*coeff,prod(72+2*j for j in range(m))*multi)-sum(val(c) for c in roots)
 else:
  assert tag[0]=='marg';i,v=tag[1:];val=lambda c:int(c[i]==v)
  b=F(M[abs(v)]-sum(val(c) for c in roots))
 row=[sum(val(cells[j]) for j in group) for group in groups]
 den=b.denominator;row=[v*den for v in row];b=int(b*den)
 g=gcd(*row,b) or 1;row=[v//g for v in row];b//=g
 rows.append(row);rhs.append(b)
matrix=sp.Matrix(rows);_,pivots=matrix.T.rref();ids=list(pivots)
print('EXACT RANK',len(ids),'rows',ids,flush=True)
A=np.array([rows[i] for i in ids],dtype=float);b=np.array([rhs[i] for i in ids],dtype=float)
sc=np.max(np.abs(A),axis=1);objective=[sum(cells[j] in map(tuple,offsets) for j in gg) for gg in groups]
res=linprog(objective,A_eq=A/sc[:,None],b_eq=b/sc/N,bounds=(0,None),method='highs')
assert res.success,res.message
ys=[F(float(y/scale)).limit_denominator(10**9) for y,scale in zip(res.eqlin.marginals,sc)]
slack=[F(v)-sum(ys[i]*rows[ii][j] for i,ii in enumerate(ids)) for j,v in enumerate(objective)]
value=sum(y*rhs[ii] for y,ii in zip(ys,ids))
print('DUAL approx',list(map(str,ys)),flush=True)
print('slack min',min(slack),'value',value,'float',float(value),flush=True)
if min(slack)<0:
 # High precision rounding plus an exact correction on the degree-zero row.
 D=10**24
 ys=[F(round(float(y/scale)*D),D) for y,scale in zip(res.eqlin.marginals,sc)]
 slack=[F(v)-sum(ys[i]*rows[ii][j] for i,ii in enumerate(ids)) for j,v in enumerate(objective)]
 z=ids.index(T['tags'].index(['mom',[0,0,0,0]]))
 delta=max([F(0)]+[-s/F(rows[ids[z]][j]) for j,s in enumerate(slack)])
 ys[z]-=delta
 slack=[F(v)-sum(ys[i]*rows[ii][j] for i,ii in enumerate(ids)) for j,v in enumerate(objective)]
 value=sum(y*rhs[ii] for y,ii in zip(ys,ids))
assert min(slack)>=0
bound=-((-value.numerator)//value.denominator)
assert bound>=743522304,(bound,value)
certificate={
 'dimension':69,'population':N,'norm_squared':8,'gram':W['gram'],
 'offsets':offsets,'heights':list(map(str,heights)),'alpha_squared':list(map(str,alpha2)),
 'cells':W['cells'],'groups':groups,
 'tags':[T['tags'][i] for i in ids],
 'rows':[[[j,v] for j,v in enumerate(rows[i]) if v] for i in ids],
 'rhs':[rhs[i] for i in ids],
 'multipliers':list(map(str,ys)),'objective':objective,
 'exact_lower_bound':str(value),'bound':bound,'slack':list(map(str,slack)),
 'geometry_gaps':gaps,
 'comparison':689452922,'gain':bound-689452922,
 'method':'D4 orientation-averaged profile populations; rational one-coordinate heads and scaled projected bodies',
}
(p/'head69_certificate.json').write_text(json.dumps(certificate,indent=2))
print('EXACT CERTIFICATE',bound,'value',value,'gain',bound-689452922,flush=True)
