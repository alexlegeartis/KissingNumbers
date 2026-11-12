"""New P48 normal-lattice models, numerical discovery only.

BCC3 has dual profiles h in D3 and p=h/(2 sqrt(2)).  Orthogonal frames
have integral profiles h and p=h/sqrt(6). Signed-permutation orientation
averaging is legitimate without assuming ambient-lattice automorphisms.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from itertools import product
from math import prod
from fractions import Fraction as F
import json,sys,subprocess
import numpy as np
from scipy.optimize import linprog

HERE=Path(__file__).parent
N=52416000
def df(n):return prod(range(n,0,-2))

def build(name,k):
    if name=='bcc':
        assert k==3
        cells=[v for v in product(range(-4,5),repeat=k) if sum(v)%2==0 and sum(map(abs,v))<=6]
        roots=list(product((-4,4),repeat=k)); den=8; moment_scale=48
    elif name=='orth':
        cells=list(product(range(-3,4),repeat=k))
        roots=[]
        for i in range(k):
            for s in (-6,6):v=[0]*k;v[i]=s;roots.append(v)
        den=6;moment_scale=36
    else:raise ValueError(name)
    orbits={}
    for i,v in enumerate(cells):orbits.setdefault(tuple(sorted(map(abs,v))),[]).append(i)
    exps=[e for e in product(range(6),repeat=k) if sum(e)<=5 and list(e)==sorted(e)]
    rows=[];rhs=[]
    for e in exps:
        s=sum(e)
        M=F(moment_scale**s*prod(df(2*a-1) for a in e),prod(48+2*j for j in range(s)))
        rows.append([float(F(sum(prod(cells[i][j]**(2*e[j]) for j in range(k)) for i in ids),4**(2*s))) for ids in orbits.values()])
        excluded=F(sum(prod(v[j]**(2*e[j]) for j in range(k)) for v in roots),N)
        rhs.append(float((M-excluded)/4**(2*s)))
    A=np.array(rows);b=np.array(rhs)
    populations=[]
    for i in range(len(orbits)):
        c=np.zeros(len(orbits));c[i]=1
        for sign in (-1,1):
            sol=linprog(sign*c,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
            if not sol.success:raise RuntimeError(sol.message)
            populations.append(sol.x*N)
    mean=np.mean(populations,axis=0)
    groups=np.zeros(len(cells),dtype=int)
    for i,ids in enumerate(orbits.values()):groups[ids]=i
    cells=np.array(cells,dtype=int)
    z=np.sqrt(3/den)*cells;q=np.sum(z*z,axis=1)
    print('POPULATIONS',name,k,'profiles',len(cells),'orbits',len(orbits),'moments',len(exps),flush=True)
    print([(key,len(ids),round(mean[i],3)) for i,(key,ids) in enumerate(orbits.items())],flush=True)
    return dict(cells=cells,groups=groups,z=z,q=q,w=mean[groups],A=A,b=b,orbits=orbits)

def run(name,k,den=16,qmax=7.5,seconds=45):
    S=build(name,k);vertices=[]
    for i in range(len(S['cells'])):
        q=S['q'][i]
        if q>qmax+1e-9 or q>=18-1e-9 or S['w'][i]<1:continue
        r=np.sqrt(q)
        heights=sorted(set([j/den for j in range(-int(r*den+1e-9),int(r*den+1e-9)+1)]+[-r,r]))
        for h in heights:
            alpha=np.sqrt((18-h*h)/(18-q))
            vertices.append((i,alpha,h,S['w'][i],*S['z'][i]))
    stem=f'{name}{k}_d{den}_q{qmax:g}'
    source=HERE/f'{stem}.txt'
    with source.open('w') as f:
        f.write(f'{len(vertices)} {k}\n')
        for row in vertices:f.write(' '.join(map(str,row))+'\n')
    print('VERTICES',len(vertices),flush=True)
    subprocess.run([str(HERE/'head_search'),str(source),str(HERE/f'{stem}_best.json'),'200000',str(seconds)],check=True)
    r=json.loads((HERE/f'{stem}_best.json').read_text());c=np.zeros(S['A'].shape[1])
    for j in r['vertices']:c[S['groups'][vertices[j][0]]]+=1
    sol=linprog(c,A_eq=S['A'],b_eq=S['b'],bounds=(0,None),method='highs')
    assert sol.success
    r.update(model=name,k=k,dimension=49-k,numerical_lower_bound=sol.fun*N,
      slices=[dict(profile=S['cells'][vertices[j][0]].tolist(),head=vertices[j][2],alpha=vertices[j][1],mean_population=vertices[j][3]) for j in r['vertices']],
      status='Numerical discovery only; exact embedding, geometry and dual verification are required')
    (HERE/f'{stem}_result.json').write_text(json.dumps(r,indent=2))
    print(json.dumps(r,indent=2),flush=True)
    return r

if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]) if len(sys.argv)>3 else 16,float(sys.argv[4]) if len(sys.argv)>4 else 7.5,float(sys.argv[5]) if len(sys.argv)>5 else 45)
