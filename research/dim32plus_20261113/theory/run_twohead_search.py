"""D4 + two head-coordinate discovery in dimension70."""
from head_system import build,bound,N
from pathlib import Path
from scipy.optimize import linprog
import numpy as np,json,sys,subprocess
k=4;r=2;folder=Path(__file__).parent
tag='_refined' if len(sys.argv)>1 else ''
S=build(k);A=S['A'];b=S['b'];population=[]
for i in range(A.shape[1]):
    e=np.zeros(A.shape[1]);e[i]=1
    for sign in [-1,1]:
        population.append(linprog(sign*e,A_eq=A,b_eq=b,bounds=(0,None),method='highs').x*N)
weight=np.mean(population,axis=0)
H=np.linalg.inv(S['gram']);z=S['cells']@np.linalg.cholesky(H);q=np.sum(z*z,axis=1)
gammas=([.125,.25,.375,.5,.625,.75,.875,1.,np.sqrt(1/3),np.sqrt(.5),np.sqrt(2/3),np.sqrt(.75),2-np.sqrt(3)] if tag else [.25,.5,.75,.875,1.])
angles=12 if tag else 24
directions=[(0.,0.)]+[(gamma*np.cos(t),gamma*np.sin(t)) for gamma in gammas for t in np.arange(angles)*2*np.pi/angles]
vertices=[]
for i,c in enumerate(S['cells']):
    if q[i]>1.01 or weight[S['groups'][i]]<1 or (tag and q[i]<1e-10):continue
    for direction in ([(0.,0.)] if q[i]<1e-10 else directions):
        h=np.array(direction)*np.sqrt(q[i]);alpha=np.sqrt((8-h@h)/(8-q[i]));vertices.append((i,alpha,weight[S['groups'][i]],*h,*z[i]))
source=folder/f'twoheads_k4{tag}.txt'
with source.open('w') as f:
    f.write(f'{len(vertices)} {k} {r}\n')
    for v in vertices:f.write(' '.join(map(str,v))+'\n')
print('vertices',len(vertices),flush=True)
subprocess.run([str(folder/'multihead_search'),str(source),str(folder/f'twoheads_k4{tag}_best.json'),'50000'],check=True)
result=json.loads((folder/f'twoheads_k4{tag}_best.json').read_text());chosen=[vertices[j][0] for j in result['vertices']]
result['k']=k;result['dimension']=72-k+r;result['numerical_lower_bound']=bound(S,chosen)
result['slices']=[{'profile':S['cells'][i].tolist(),'head':list(vertices[j][3:5]),'alpha':vertices[j][1],'weight':vertices[j][2]} for j,i in zip(result['vertices'],chosen)]
result['status']='NUMERICAL DISCOVERY ONLY; geometry and population require exact certification'
(folder/f'twoheads_k4{tag}_result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
