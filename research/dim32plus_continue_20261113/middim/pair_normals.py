"""Adaptive 47D compression with two norm-six normals of arbitrary overlap.

Populations are obtained from all bivariate moments through degree eleven.
The remaining one-dimensional population interval is treated adversarially.
Numerical search output is not an exact certificate.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod,comb,factorial
import numpy as np
import sympy as S
import json,sys,subprocess
HERE=Path(__file__).parent
N=52416000
def df(n):return prod(range(n,0,-2))
def moment(a,b,c):
    if (a+b)%2:return F(0)
    k=(a+b)//2
    z=sum(comb(a,t)*comb(b,t)*factorial(t)*df(a-t-1)*df(b-t-1)*6**((a+b-2*t)//2)*c**t for t in range(min(a,b)+1) if (a-t)%2==0 and (b-t)%2==0)
    return F(6**k*z,prod(48+2*i for i in range(k)))
def run(c,den=48,seconds=60):
    assert c in (-2,-1,0)
    cells=list(product(range(-3,4),repeat=2));exps=[(a,b) for a,b in product(range(7),repeat=2) if a+b<=11]
    roots=[(6,c),(-6,-c),(c,6),(-c,-6)]
    A=S.Matrix([[x**a*y**b for x,y in cells] for a,b in exps])
    b=S.Matrix([S.Rational(N)*S.Rational(moment(a,b,c))-sum(x**a*y**b for x,y in roots) for a,b in exps])
    counts,t=A.gauss_jordan_solve(b);assert len(t)==1;t=t[0]
    lo=S.Rational(0);hi=S.oo
    for f in counts:
        slope=f.diff(t);offset=f.subs(t,0)
        if slope>0:lo=max(lo,-offset/slope)
        elif slope<0:hi=min(hi,-offset/slope)
        else:assert offset>=0
    lo=S.ceiling(lo);hi=S.floor(hi);assert lo<=hi
    print('MOMENTS overlap',c,'parameter interval',lo,hi,flush=True)
    weights=np.array([float(f.subs(t,(lo+hi)/2)) for f in counts]);p=np.array(cells,dtype=int)
    G=np.array([[6,c],[c,6]]);z=p@np.linalg.cholesky(3*np.linalg.inv(G));q=np.sum(z*z,axis=1)
    vertices=[]
    for i in range(len(cells)):
        if q[i]>=18-1e-8 or weights[i]<1:continue
        r=np.sqrt(q[i]);heights=sorted(set([j/den for j in range(-int(r*den+1e-8),int(r*den+1e-8)+1)]+[-r,r]))
        for h in heights:vertices.append((i,np.sqrt((18-h*h)/(18-q[i])),h,weights[i],*z[i]))
    stem=f'pair{c}_d{den}';source=HERE/f'{stem}.txt'
    with source.open('w') as f:
        f.write(f'{len(vertices)} 2\n')
        for row in vertices:f.write(' '.join(map(str,row))+'\n')
    subprocess.run([str(HERE/'head_search'),str(source),str(HERE/f'{stem}_best.json'),'200000',str(seconds)],check=True)
    r=json.loads((HERE/f'{stem}_best.json').read_text());objective=sum(counts[vertices[j][0]] for j in r['vertices'])
    r.update(overlap=c,dimension=47,numerical_lower_bound=int(min(objective.subs(t,lo),objective.subs(t,hi))),
      parameter_interval=[int(lo),int(hi)],population_expression=str(objective),
      slices=[dict(profile=cells[vertices[j][0]],head=vertices[j][2],alpha=vertices[j][1],population=str(counts[vertices[j][0]])) for j in r['vertices']],
      status='Numerical discovery only; exact geometry and embedding required')
    (HERE/f'{stem}_result.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2),flush=True)
if __name__=='__main__':run(int(sys.argv[1]),int(sys.argv[2]) if len(sys.argv)>2 else 48,float(sys.argv[3]) if len(sys.argv)>3 else 60)
