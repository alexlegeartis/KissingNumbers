"""One-head A4 compression search; does not update any claimed bound."""
from pathlib import Path
from a4_model import build,solve,N
import numpy as np,json,subprocess
from scipy.optimize import linprog
P=Path(__file__).parent;M=build();population=[]
for j in range(len(M['groups'])):
 a=np.zeros(len(M['groups']));a[j]=1
 for sign in [-1,1]:
  r=linprog(sign*a*M['colscale'],A_eq=M['As'],b_eq=M['bs'],bounds=(0,None),method='highs',options={'threads':1});assert r.success
  population.append(r.x*M['colscale']*N)
w=np.mean(population,axis=0)
gammas=sorted(set([0]+[s*x for s in [-1,1] for x in [.125,.25,.375,.5,.625,.75,.875,1,np.sqrt(1/3),np.sqrt(.5),np.sqrt(2/3),np.sqrt(.75),2-np.sqrt(3)]]))
verts=[]
for i,c in enumerate(M['cells']):
 q=M['q'][i]
 if q>1.5+1e-8 or w[M['groupid'][i]]<1:continue
 for g in ([0] if q<1e-8 else gammas):
  h=g*np.sqrt(q);alpha=np.sqrt((8-h*h)/(8-q));verts.append((i,alpha,h,w[M['groupid'][i]],*M['z'][i]))
with (P/'a4_vertices.txt').open('w') as f:
 f.write(f'{len(verts)} 5\n')
 for row in verts:f.write(' '.join(map(str,row))+'\n')
print('vertices',len(verts),flush=True)
engine=Mpath=P.parents[0]/'dim32plus_20261113/theory/head_search'
engine=P.resolve().parents[2]/'research/dim32plus_20261113/theory/head_search'
subprocess.run([str(engine),str(P/'a4_vertices.txt'),str(P/'a4_best.json'),'1000000'],check=True)
r=json.loads((P/'a4_best.json').read_text());selected=[verts[j][0] for j in r['vertices']]
value,lp,obj=solve(M,selected)
r.update(numerical_population_lower=value,objective=obj.tolist(),slices=[{'profile':M['cells'][i].tolist(),'normal_projection':M['z'][i].tolist(),'q':M['q'][i],'height':verts[j][2],'alpha':verts[j][1]} for j,i in zip(r['vertices'],selected)],status='NUMERICAL RESEARCH ONLY')
(P/'a4_result.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
