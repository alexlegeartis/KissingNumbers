"""Generate a discretized head-compression compatibility graph for discovery."""
from head_system import build,bound,N
from pathlib import Path
from scipy.optimize import linprog
import numpy as np,json,sys,subprocess
k=int(sys.argv[1]);folder=Path(__file__).parent
S=build(k);A=S['A'];b=S['b'];population=[]
for i in range(A.shape[1]):
    e=np.zeros(A.shape[1]);e[i]=1
    for sign in [-1,1]:
        population.append(linprog(sign*e,A_eq=A,b_eq=b,bounds=(0,None),method='highs').x*N)
weight=np.mean(population,axis=0)
H=np.linalg.inv(S['gram']);z=S['cells']@np.linalg.cholesky(H);q=np.sum(z*z,axis=1)
gammas=sorted(set([0]+[s*x for s in [-1,1] for x in [0.125,.25,.375,.5,.625,.75,.875,1.,np.sqrt(1/3),np.sqrt(.5),np.sqrt(2/3),np.sqrt(.75),2-np.sqrt(3)]]))
vertices=[]
for i,c in enumerate(S['cells']):
    if q[i]>2.01 or weight[S['groups'][i]]<1:continue
    for gamma in ([0] if q[i]<1e-10 else gammas):
        h=gamma*np.sqrt(q[i]);alpha=np.sqrt((8-h*h)/(8-q[i]));vertices.append((i,alpha,h,weight[S['groups'][i]],*z[i]))
source=folder/f'heads_k{k}.txt'
with source.open('w') as f:
    f.write(f'{len(vertices)} {k}\n')
    for v in vertices:f.write(' '.join(map(str,v))+'\n')
print('vertices',len(vertices),flush=True)
subprocess.run([str(folder/'head_search'),str(source),str(folder/f'heads_k{k}_best.json'),'30000'],check=True)
r=json.loads((folder/f'heads_k{k}_best.json').read_text());chosen=[vertices[j][0] for j in r['vertices']]
r['k']=k;r['dimension']=73-k;r['numerical_lower_bound']=bound(S,chosen);r['slices']=[{'profile':S['cells'][i].tolist(),'h':vertices[j][2],'alpha':vertices[j][1],'weight':vertices[j][3]} for j,i in zip(r['vertices'],chosen)]
r['status']='NUMERICAL DISCOVERY ONLY; geometry and population require exact certification'
(folder/f'heads_k{k}_result.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
