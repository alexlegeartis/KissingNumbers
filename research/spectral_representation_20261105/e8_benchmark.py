"""Recover the known 240-point E8 code from its supplied contact graph."""
from itertools import combinations,product
from pathlib import Path
import sys,time,numpy as np
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root))
from discovery.core import RunContext
from discovery.engines.spectral_graph_search import search
from discovery.certification import verify_artifact
X=[]
for i,j in combinations(range(8),2):
 for a,b in product((-2,2),repeat=2):
  v=np.zeros(8,np.int64);v[i]=a;v[j]=b;X.append(v)
for signs in product((-1,1),repeat=8):
 if sum(x<0 for x in signs)%2==0:X.append(np.array(signs,np.int64))
X=np.array(X);assert X.shape==(240,8) and np.all((X*X).sum(1)==8);A=(X@X.T==4).astype(int);np.fill_diagonal(A,0)
ctx=RunContext(root,{'max_rank_work':500000000,'solver_seconds':30});t=time.time();r=search([],{'adjacency':A.tolist(),'max_dimension':8,'max_sectors':1,'include_antipodes':False,'seconds':30},ctx);assert r.facts['cardinality']==240 and r.facts['dimension']==8 and verify_artifact(r,ctx).status=='exact';print('PASS 240 points in dimension 8 from the known E8 contact graph; seconds',round(time.time()-t,1))
