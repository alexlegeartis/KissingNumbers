import numpy as np,time,os
from scipy.sparse import coo_matrix
from scipy.optimize import milp,Bounds,LinearConstraint

def optimize(C,blockers,norm=128,seconds=90):
 st=time.monotonic();n=len(C);old=np.unique([j for b in blockers for j in b]);mapping={int(x):n+i for i,x in enumerate(old)}
 r=[];c=[];v=[];upper=[]
 for i,bs in enumerate(blockers):
  for j in bs:
   k=len(upper);r.extend([k,k]);c.extend([i,mapping[j]]);v.extend([1,-1]);upper.append(0)
 C=C.astype(np.float32);candedges=0
 for at in range(0,n,256):
  M=C[at:at+256]@C.T
  for loc in range(len(M)):
   i=at+loc;js=np.flatnonzero(M[loc,i+1:]>norm//2)+i+1
   for j in js:
    k=len(upper);r.extend([k,k]);c.extend([i,int(j)]);v.extend([1,1]);upper.append(1)
   candedges+=len(js)
 print('MODEL',n,len(old),candedges,'seconds',time.monotonic()-st,flush=True)
 mat=coo_matrix((v,(r,c)),shape=(len(upper),n+len(old))).tocsc()
 ans=milp(np.r_[-np.ones(n),np.ones(len(old))],integrality=np.ones(n+len(old)),bounds=Bounds(0,1),constraints=LinearConstraint(mat,-np.inf,upper),options={'time_limit':seconds,'mip_rel_gap':0,'threads':int(os.environ.get('RESEARCH_THREADS','3'))})
 add=[] if ans.x is None else np.flatnonzero(ans.x[:n]>.5).tolist();remove=[] if ans.x is None else old[np.flatnonzero(ans.x[n:]>.5)].tolist()
 return {'gain':len(add)-len(remove),'added':add,'removed':remove,'status':int(ans.status),'message':str(ans.message),'new_conflicts':candedges,'numerical_gain_upper_bound':-float(ans.mip_dual_bound) if getattr(ans,'mip_dual_bound',None) is not None else None,'scope':'Numerical optimizer, exact independent replay required for positive result.'}
