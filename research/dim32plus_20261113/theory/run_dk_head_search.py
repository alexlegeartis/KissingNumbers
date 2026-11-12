"""One-head discovery using root's independent orthogonal D_k model."""
from pathlib import Path
import sys,json,subprocess
import numpy as np
from scipy.optimize import linprog
folder=Path(__file__).parent
sys.path.insert(0,str(folder.parent/'gamma'))
from dk_affine import build,solve_obj,N,serialize_model
k=int(sys.argv[1]);qmax=float(sys.argv[2]) if len(sys.argv)>2 else 1.5
M=build(k);population=[]
for i in range(len(M['groups'])):
    e=np.zeros(len(M['groups']));e[i]=1
    for sign in [-1,1]:
        r=linprog(sign*e*M['colscale'],A_eq=M['As'],b_eq=M['bs'],bounds=(0,None),method='highs',options={'threads':1})
        assert r.success
        population.append(r.x*M['colscale']*N)
weight=np.mean(population,axis=0)
cells=np.array(sorted(M['gid']),dtype=int);z=cells/4;q=np.sum(z*z,axis=1)
gammas=sorted(set([0]+[s*x for s in [-1,1] for x in [.125,.25,.375,.5,.625,.75,.875,1.,np.sqrt(1/3),np.sqrt(.5),np.sqrt(2/3),np.sqrt(.75),2-np.sqrt(3)]]))
vertices=[]
for i,c in enumerate(cells):
    if q[i]>qmax+.001 or weight[M['gid'][tuple(c)]]<1:continue
    for gamma in ([0] if q[i]<1e-10 else gammas):
        h=gamma*np.sqrt(q[i]);alpha=np.sqrt((8-h*h)/(8-q[i]));vertices.append((i,alpha,h,weight[M['gid'][tuple(c)]],*z[i]))
source=folder/f'dheads_k{k}.txt'
with source.open('w') as f:
    f.write(f'{len(vertices)} {k}\n')
    for v in vertices:f.write(' '.join(map(str,v))+'\n')
print('vertices',len(vertices),flush=True)
subprocess.run([str(folder/'head_search'),str(source),str(folder/f'dheads_k{k}_best.json'),'50000'],check=True)
r=json.loads((folder/f'dheads_k{k}_best.json').read_text());chosen=[vertices[j][0] for j in r['vertices']]
r['k']=k;r['dimension']=73-k;r['numerical_lower_bound']=solve_obj(M,cells[chosen])[0];r['slices']=[{'profile_h':cells[i].tolist(),'h':vertices[j][2],'alpha':vertices[j][1],'weight':vertices[j][3]} for j,i in zip(r['vertices'],chosen)]
r['status']='NUMERICAL DISCOVERY ONLY; geometry and population require exact certification'
(folder/f'dheads_k{k}_result.json').write_text(json.dumps(r,indent=2));(folder/f'dheads_k{k}_model.json').write_text(json.dumps(serialize_model(M)));print(json.dumps(r,indent=2))
