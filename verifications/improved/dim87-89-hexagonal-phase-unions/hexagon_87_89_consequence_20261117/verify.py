import numpy as np,json,itertools
from pathlib import Path
P=Path(__file__).resolve().parent
out=[]
for m,M in [(3,9),(4,12),(5,16)]:
 roots=[]
 for i,j in itertools.combinations(range(m),2):
  for a,b in itertools.product([-1,1],repeat=2):
   x=np.zeros(m,dtype=np.int64);x[i]=a;x[j]=b;roots.append(x)
 R=np.array(roots);D=R@R.T
 assert len(R)==2*m*(m-1) and len(np.unique(R,axis=0))==len(R)
 assert np.all(np.diag(D)==2)and D[np.triu_indices(len(R),1)].max()==1
 # Scaled norm45 root v=sqrt(45/2)*r versus old norm3 head: Cauchy<=sqrt135<45/2.
 assert 4*135<45**2
 dense=2*M*3**18;sparse=37468*7812;total=dense+sparse+len(R)
 out.append(dict(dimension=84+m,dense=dense,sparse=sparse,pure_head=len(R),total=total))
assert [x['total'] for x in out]==[7266268830,9590791776,12690155704]
(P/'certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
