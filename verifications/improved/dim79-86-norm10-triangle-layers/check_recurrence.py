from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,hashlib
import sympy as s
import sys
n=int(sys.argv[1]);P=Path(__file__).resolve().parent;S=P;x=json.load(open(S/f'exact{n}_marginal.json'));y=json.load(open(S/f'blocked{n}_marginal.json'));D=x['degree'];T=Q(x['triangle_upper']);assert T.denominator==1 and T>=0
# Derive transformed recurrence, without assuming producer coefficient formula.
d,t,g=s.symbols('d t g');den=d*d+1-2*t
trans=s.cancel((d+1)*(1+(d*d-d-2*t)*(g+1)/d)/den-1)
assert s.simplify(s.diff(trans,g)-(1-2*(d+t)/(d*den)))==0
assert s.simplify(trans.subs(g,0)-(d*d-d-2*t)/(d*den))==0
# Ratios used for monotonicity decrease above the positive-denominator threshold.
r=(d+1)/den
assert s.simplify(s.diff(r,d)+(d*d+2*d-1+2*t)/den**2)==0
# (1+t/d)/den decreases since both nonnegative factors do.
start=y['blocks'][0]['first_degree'];last=start-1;lower=Q(0);checks=[]
for block in y['blocks']:
 a=block['first_degree'];b=block['last_degree'];assert a==last+1 and b>=a and b<=D
 assert a*(a-1)>=2*T and a>0
 den=a*a+1-2*T;assert den>0
 B=(a*a-a-2*T)/(b*den);C=2*(a+T)/(a*den)
 assert Q(block['B'])==B and Q(block['C'])==C and B>=0 and 0<C<1
 L=b-a+1;lower=max(Q(0),1-L*C)*lower+B*L/(1+L*C)
 checks.append(dict(a=a,b=b,length=L));last=b
assert len(checks)==40 and last==D
N=15281788354560//2*(16*comb(n-72,4));assert x['tail_count']==16*comb(n-72,4) and N==y['graph_vertices']
assert lower+1==Q(y['factor']);alpha=N*(lower+1)/(D+1);assert alpha==Q(y['expected_pairs'])
pairs=-(-alpha.numerator//alpha.denominator);assert pairs==y['pairs']
bound=6218175600+2*pairs+2*(n-72)*(n-73);assert bound==y['bound']
report=dict(status='PASS transformed recurrence and every rational block; conditional on graph/theorem premises audited separately',dimension=n,blocks=40,degree=D,triangles=str(T),graph_vertices=N,pairs=pairs,bound=bound,block_start=start,source_sha256={f:hashlib.sha256((S/f).read_bytes()).hexdigest()for f in [f'exact{n}_marginal.json',f'blocked{n}_marginal.json']})
(P/f'RESULT{n}_marginal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
