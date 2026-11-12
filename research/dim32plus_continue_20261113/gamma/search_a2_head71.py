"""Bounded A2-to-one-head search using D7-averaged exact population models."""
from pathlib import Path
from fractions import Fraction as F
import json
import math
import subprocess
import numpy as np
from scipy.optimize import linprog
import extension_model as m

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
base=json.loads((HERE/'d7_compact70.json').read_text())
reps=[tuple(h) for h in base['representatives']]
sizes=np.asarray(base['orbit_sizes'],float)
A=np.asarray(base['rows'],float)/sizes[None,:]
b=np.asarray(base['rhs'],float)/m.N
scales=np.maximum(abs(b),1e-15)
options=dict(threads=1,presolve=False,primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9)
cells=[(a,b) for a in range(-4,5) for b in range(-4,5) if abs(a-b)<=4]
objectives=[];known=[]
for c in cells:
    marginal=m.a2_marginal_weights([c])
    objectives.append([float(m.objective(h,marginal)) for h in reps])
    known.append(float(sum((marginal.get(tuple(sorted(map(abs,h[:3]))),F(0)) for h in m.normal_roots(7)),F(0))))
C=np.asarray(objectives)/sizes[None,:]
known=np.asarray(known)
samples=[]
for row in C:
    for sign in (-1,1):
        lp=linprog(sign*row,A_eq=A/scales[:,None],b_eq=b/scales,bounds=(0,None),method='highs',options=options)
        assert lp.success,lp.message
        samples.append(lp.x)
weights=C@np.mean(samples,axis=0)*m.N+known
L=np.linalg.cholesky(np.array([[1/6,-1/12],[-1/12,1/6]]))
vertices=[];metadata=[]
for i,c in enumerate(cells):
    n=c[0]**2-c[0]*c[1]+c[1]**2
    q=F(n,6)
    limit=math.isqrt(n*144)
    z=np.asarray(c)@L
    for j in range(-limit,limit+1):
        u=F(j,12)
        t2=u*u/6
        assert t2<=q
        alpha=math.sqrt(float((8-t2)/(8-q)))
        vertices.append((i,alpha,float(u)/math.sqrt(6),weights[i],*z))
        metadata.append((i,u))
engine=REPO/'research/dim32plus_20261113/theory/head_search'
print('A2_HEAD71 vertices',len(vertices),'cells',len(cells),flush=True)
best=0
for attempt in range(3):
    inputfile=HERE/f'a2_head71_vertices_{attempt}.txt'
    with inputfile.open('w') as f:
        f.write(f'{len(vertices)} 2\n')
        for i,row in enumerate(vertices):
            row=list(row);row[3]=weights[row[0]]
            f.write(' '.join(map(str,row))+'\n')
    resultfile=HERE/f'a2_head71_search_{attempt}.json'
    subprocess.run([str(engine),str(inputfile),str(resultfile),'1000000'],check=True)
    raw=json.loads(resultfile.read_text())
    selected=[metadata[j] for j in raw['vertices']]
    profiles=[cells[i] for i,u in selected]
    heights=[str(u) for i,u in selected]
    inner=lambda c,d:F(2*c[0]*d[0]-c[0]*d[1]-c[1]*d[0]+2*c[1]*d[1],12)
    q=[inner(c,c) for c in profiles]
    u=list(map(F,heights));alpha2=[(8-t*t/6)/(8-r) for t,r in zip(u,q)]
    for i in range(len(profiles)):
        for j in range(i+1):
            dd=4-inner(profiles[i],profiles[j]);rr=4-u[i]*u[j]/6
            assert dd>=0 and rr>=0 and rr*rr-alpha2[i]*alpha2[j]*dd*dd>=0,(i,j)
    result=m.solve(7,71,profiles,heights,F(1,6),f'a2_head71_candidate_{attempt}.json')
    value=result.get('with_two_poles',result.get('exact_bound',result['value']))
    print('A2_HEAD71 attempt',attempt,'slices',len(selected),'value',value,'prior',2603658750,flush=True)
    if value>best:
        best=value
        (HERE/'a2_head71_best.json').write_text(json.dumps(result,indent=2)+'\n')
    # Use the exact-model worst-case population for the next geometric search.
    mass=np.asarray(result['primal'])*sizes
    weights=C@mass*m.N+known
