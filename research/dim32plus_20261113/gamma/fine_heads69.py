"""Finer one-head grid; numerical discovery, not a proof."""
import sys,json,subprocess
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
P=Path(__file__).parent
sys.path.insert(0,str(P.parent/'theory'))
from head_system import build,bound,N
S=build(4);A=S['A'];b=S['b'];pop=[]
for j in range(A.shape[1]):
    e=np.zeros(A.shape[1]);e[j]=1
    for sign in [-1,1]:pop.append(linprog(sign*e,A_eq=A,b_eq=b,bounds=(0,None),method='highs').x*N)
w=np.mean(pop,axis=0)
H=np.linalg.inv(S['gram']);z=S['cells']@np.linalg.cholesky(H);q=np.sum(z*z,axis=1)
vertices=[]
for i,c in enumerate(S['cells']):
    if q[i]>.75000001 or w[S['groups'][i]]<1:continue
    if q[i]<1e-10:heights=[0.]
    else:
        step=128 if q[i]<.51 else 32
        end=int(np.sqrt(q[i])*step+1e-9)
        heights=sorted(set([j/step for j in range(-end,end+1)]+[-np.sqrt(q[i]),np.sqrt(q[i])]))
    for h in heights:
        alpha=np.sqrt((8-h*h)/(8-q[i]));vertices.append((i,alpha,h,w[S['groups'][i]],*z[i]))
source=P/'fine_heads69.txt'
with source.open('w') as f:
    f.write(f'{len(vertices)} 4\n')
    for v in vertices:f.write(' '.join(map(str,v))+'\n')
print('VERTICES',len(vertices),flush=True)
subprocess.run([str(P.parent/'theory/head_search'),str(source),str(P/'fine_heads69_best.json'),'100000'],check=True)
r=json.loads((P/'fine_heads69_best.json').read_text());chosen=[vertices[j][0] for j in r['vertices']]
r['numerical_lower_bound']=bound(S,chosen)
r['slices']=[{'profile':S['cells'][vertices[j][0]].tolist(),'h':vertices[j][2],'alpha':vertices[j][1]} for j in r['vertices']]
r['status']='Numerical discovery only; exact verification required'
(P/'fine_heads69_result.json').write_text(json.dumps(r,indent=2))
print(json.dumps(r,indent=2))
